% This script simulates alpha of varying cycles and amplitudes superimposed
% on pink noise. Subsequently, the STANDARD BOSC procedure is run
% Resulting abundance and signal detection measures are saved 
% to compare the measures' sensitivity and specificity.

%% set up paths

    pn.root = '/PATH-TO-eBOSC-TOOLBOX/eBOSC/'; addpath(genpath(pn.root));
    pn.backgroundData  = ''; % INDICATE LOCATION OF SIMULATED BACKGROUND
    pn.out = ''; % CHOOSE OUTPUT DIRECTORY

%%  parameters

    % BOSC parameters   
    bosc.F          = 2.^[1:.125:5.25];  % setup (Whitten et al., 2011), but higher frequency resolution
    bosc.wavenumber = 6;    
    bosc.fsample    = 250;
    bosc.WLpadding   = 500; % padding to avoid edge artifacts due to WL [SPs]
    bosc.detectedPad = 250; % 'shoulder' for BOSC detected matrix to account for duration threshold
    bosc.trialPad    = 750; % complete padding (WL + shoulder)
    bosc.BGpad       = 750; % padding of segments for BG (only avoiding edge artifacts)
    
    % BOSC parameters
    numcycles  = [3];
    percentile = [.95];
    
    % alpha amplitudes
    amplitude = [0 2 4 6 8 12 16 24];
    
    % alpha cycles 
    cycles = [2 4 8 16 32 64 128 200]; % total of 14 seconds
    
    % time vector
    timeDur = 20;
    time = [.004:.004:timeDur];
    
    % set amount of repetitions
    
    repetitions = 500;
        
%%  load background

    load([pn.backgroundData 'background.mat'],'bckgrnd_filt')
    
%%  loop conditions

SignalDetection.Hits= NaN(numel(amplitude), numel(cycles), repetitions,1);
SignalDetection.Misses= NaN(numel(amplitude), numel(cycles), repetitions,1);
SignalDetection.CAs= NaN(numel(amplitude), numel(cycles), repetitions,1);
SignalDetection.FAs= NaN(numel(amplitude), numel(cycles), repetitions,1);

abundance_spec= NaN(numel(amplitude), numel(cycles), repetitions, numel(bosc.F));

for n = 1 % numcycles        
for p = 1 % percentile (actually, we could just use the .95 one)

    bosc.numcycles  = numcycles(n);
    bosc.percentile = percentile(p);

for a = 1:length(amplitude)
for c = 1:length(cycles)
        
for k = 1:repetitions
    
    %  generate alpha in the middle of the segment

    alphaCycles = cycles(c);
    alphaFreq = 10;
    alphaTime(c) = round((alphaCycles/alphaFreq),3);
    % make alpha symmetrical around the middle
    timeNew = round(alphaTime(c)/0.004,0);
    if mod(timeNew,2) ~= 0
        timeNew = timeNew + 1;
        alphaTime(c) = timeNew.*0.004;
    else alphaTime(c) = timeNew.*0.004;
    end
    timeAlpha = [.004:.004:alphaTime(c)];
    AlphaPlace = (numel(time)/2)-(numel(timeAlpha)/2)+1:(numel(time)/2)+(numel(timeAlpha)/2);
    
    % filter entire signal between 8 and 12 Hz (6th order butterworth) (not locally on alpha)
    [tmp_b,tmp_a] = butter(6, [8, 12]/(250/2), 'bandpass');
    tmp_bpsignal = filter(tmp_b,tmp_a,bckgrnd_filt(k,:));
    VarBG = var(tmp_bpsignal);
    
    targetPower = amplitude(a)*VarBG;
    amplitudeFromRMS = (sqrt(targetPower)*sqrt(2));
    alpha_sim = sin(timeAlpha*2*pi*alphaFreq)*amplitudeFromRMS;
    alpha = zeros(1,numel(time));
    alpha(1,AlphaPlace) = alpha_sim;
    
    signal = bckgrnd_filt(k,:) + alpha;
    
    %  wavelet transform
    B(k,:,:) = BOSC_tf(signal,bosc.F,bosc.fsample,bosc.wavenumber);
end

%%  BOSC - standard fit

    % background power estimation
    
    BG = [];
    for k = 1:repetitions
        BG = cat(3,BG, B(k,:,:));
    end
    BG = squeeze(BG);
    
    [pv,~] = BOSC_bgfit(bosc.F,BG);
    
    mp = 10.^(polyval(pv,log10(bosc.F))); 
    
    % thresholds
    [pt,dt] = BOSC_thresholds(bosc.fsample,bosc.percentile,bosc.numcycles,bosc.F,mp);
    
    for k = 1:repetitions
    
        display([num2str(a) '/8, ' num2str(c) '/8, ' num2str(k) '/500'])
        
        % oscillation detection
        detected = zeros(35,numel(time)-2*bosc.WLpadding);
        for f = 1:length(bosc.F)
            detected(f,:) = BOSC_detect(squeeze(B(k,f,bosc.WLpadding+1:end-bosc.WLpadding))',pt(f),dt(f),bosc.fsample);
        end; clear f

        % detected
        abundance_spec(a,c,k,:) = mean(detected(:,bosc.detectedPad+1:end-bosc.detectedPad),2);
        
        % calculate abundance within alpha range
        
        detected_meanAlpha = mean(detected(find(bosc.F >= 8 & bosc.F <= 12),bosc.detectedPad+1:end-bosc.detectedPad),1);
        detected_meanAlpha(find(detected_meanAlpha>0)) = 1; % binary vector of alpha detected or not across frequencies in alpha range
        pepisode_meanAlpha = numel(find(detected_meanAlpha==1))./numel(detected_meanAlpha);
        
        meanPepisode_Alpha = mean(mean(abundance_spec(a,c,k,find(bosc.F >= 8 & bosc.F <= 12))));
        
        % Note that the above indices are not equivalent! 
        % pepisode_meanAlpha is the pepisode of any sample points detected
        % in the alpha range; meanPepisode_Alpha is the arithmetic mean of
        % the pepisode values for each frequency in the alpha range
        
        % extract sample points with alpha detected       
       
        numberOfEffectivePoints = (numel(time)-2*bosc.trialPad);
        
        % create unique detected vector
        alpha_locs = find(detected_meanAlpha);
        % extract abundance
        abundance_PepMeanAlpha(a,c,k,1) = pepisode_meanAlpha;
        abundance_meanPep(a,c,k,1) = meanPepisode_Alpha;
    
        if a > 1
            % identify data points in final episode space that were
            % simulated as rhythmic. Note that this gives an
            % accurate amount of total signals as well.
            alpha_sim_locs = zeros(1,numel(time));
            alpha_sim_locs(AlphaPlace) = 1;
            alpha_sim_locs = find(alpha_sim_locs(bosc.trialPad+1:end-bosc.trialPad) == 1);
        else
            alpha_sim_locs = [];
        end

        % encode amount of alpha/non-alpha
        Amount.Alpha(c) = numel(alpha_sim_locs);
        Amount.NoAlpha(c) = numberOfEffectivePoints-numel(alpha_sim_locs);

        Hits = intersect(alpha_locs, alpha_sim_locs); % detected & simulated
        SignalDetection.Hits(1,a,c,k,1) = numel(Hits);
        SignalDetection.HitRate(1,a,c,k,1) = numel(Hits)/numel(alpha_sim_locs);
        FAs = setdiff(alpha_locs, alpha_sim_locs); % detected, not simulated
        SignalDetection.FAs(1,a,c,k,1) = numel(FAs);
        SignalDetection.FARate(1,a,c,k,1) = numel(FAs)/(numberOfEffectivePoints-numel(alpha_sim_locs));
        Misses = setdiff(alpha_sim_locs, alpha_locs); % simulated, not detected
        SignalDetection.Misses(1,a,c,k,1) = numel(Misses);
        SignalDetection.MissRate(1,a,c,k,1) = numel(Misses)/numel(alpha_sim_locs);
        CRs = setdiff(1:numberOfEffectivePoints, unique([alpha_sim_locs, alpha_locs])); % not simulated, not detected
        SignalDetection.CRs(1,a,c,k,1) = numel(CRs);
        SignalDetection.CRRate(1,a,c,k,1) = numel(CRs)/(numberOfEffectivePoints-numel(alpha_sim_locs));

        %% extract parameters
        
        curTF = squeeze(B(k,:,bosc.WLpadding+1:end-bosc.WLpadding));
        maskedTF = detected.*squeeze(B(k,:,bosc.WLpadding+1:end-bosc.WLpadding));
        maskedTF(maskedTF==0) = NaN;
        targetRows = find(bosc.F >= 8 & bosc.F <= 12);
        [maxVal, maxRow] = max(nanmean(maskedTF(targetRows,:),2));
        [maxValo, maxRowo] = max(nanmean(curTF(targetRows,:),2));
        
        SignalDetection.oFreq(1,a,c,k,1) = bosc.F(targetRows(1)+maxRowo-1);
        SignalDetection.oAmp(1,a,c,k,1) = sqrt(maxValo);
        
        if ~isempty(maxVal) && ~isnan(maxVal)
            SignalDetection.Freq(1,a,c,k,1) = bosc.F(targetRows(1)+maxRow-1);
            SignalDetection.Amp(1,a,c,k,1) = sqrt(maxVal);
            SignalDetection.oFreq(1,a,c,k,1) = bosc.F(targetRows(1)+maxRowo-1);
            SignalDetection.oAmp(1,a,c,k,1) = sqrt(maxValo);

            SignalDetection.fitBG(1,a,c,k,1) = sqrt(mp(targetRows(1)+maxRow-1));
            SignalDetection.PT(1,a,c,k,1) = pt(targetRows(1)+maxRow-1);
            SignalDetection.BGdiff(1,a,c,k,1) = sqrt(maxVal)-sqrt(mp(targetRows(1)+maxRow-1));
            SignalDetection.BGrel(1,a,c,k,1) = (sqrt(maxVal)-sqrt(mp(targetRows(1)+maxRow-1)))/sqrt(mp(targetRows(1)+maxRow-1));
            SignalDetection.abn(1,a,c,k,1) = pepisode_meanAlpha;
            SignalDetection.abn_THG(1,a,c,k,1) = pepisode_meanAlpha;
            SignalDetection.NaNmat(1,a,c,k,1) = 0;
        else 
            SignalDetection.Freq(1,a,c,k,1) = NaN;
            SignalDetection.Amp(1,a,c,k,1) = NaN;
            SignalDetection.oFreq(1,a,c,k,1) = NaN;
            SignalDetection.oAmp(1,a,c,k,1) = NaN;

            SignalDetection.fitBG(1,a,c,k,1) = NaN;
            SignalDetection.PT(1,a,c,k,1) = NaN;
            SignalDetection.BGdiff(1,a,c,k,1) = NaN;
            SignalDetection.BGrel(1,a,c,k,1) = NaN;
            SignalDetection.abn(1,a,c,k,1) = NaN;
            SignalDetection.abn_THG(1,a,c,k,1) = NaN;
            SignalDetection.NaNmat(1,a,c,k,1) = NaN;
        end
        
    end; clear k; % iterations
    
end; clear c; % cycles

end; clear a; % amplitudes

end; clear p; % percentile

end; clear n; % numcycles

save([pn.out, 'REDSimulation_standardBOSC_170630.mat'],'-v7.3');
