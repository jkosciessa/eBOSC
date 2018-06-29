function C_simulation_extendedBOSC(variant)

% This function simulates alpha of varying cycles and amplitudes superimposed
% on pink noise. Subsequently, the extended BOSC procedure is run using
% either the THG or FWHM post-processing. Resulting abundance and signal
% detection measures are saved for comparison between the measures'
% sensitivity and specificity.

% The input 'variant' is expected to be a string (A to H), which will call
% the detection post-processing with different parameters (see 175ff.).

%% set up paths

pn.root = '/PATH-TO-eBOSC-TOOLBOX/eBOSC/'; addpath(genpath(pn.root));
pn.backgroundData  = ''; % INDICATE LOCATION OF SIMULATED BACKGROUND
pn.out = ''; % CHOOSE OUTPUT DIRECTORY

%% encode run date

cfg.runDate = date;

%% simulation parameters

cfg.simParams.amplitude     = [0 2 4 6 8 12 16 24];                         % simulated signal power
cfg.simParams.cycles        = [2 4 8 16 32 64 128 200];                     % simulated signal durations [in cycles]; total of 14 seconds
cfg.simParams.segmentDur    = 20;
cfg.simParams.time          = [.004:.004:cfg.simParams.segmentDur];         % time vector of complete segment
cfg.simParams.repetitions   = 500;                                          % amount of repetitions

%%  eBOSC parameters

cfg.eBOSC.F                 = 2.^[1:.125:5.25];                             % setup (Whitten et al., 2011), but higher frequency resolution
cfg.eBOSC.wavenumber        = 6;
cfg.eBOSC.ncyc              = repmat(3, 1, numel(cfg.eBOSC.F));
cfg.eBOSC.percentile        = .95;
cfg.eBOSC.fsample           = 250;
cfg.eBOSC.WLpadding         = 500;                                          % padding to avoid edge artifacts due to WL [SPs]
cfg.eBOSC.detectedPad       = 250;                                          % 'shoulder' for BOSC detected matrix to account for duration threshold
cfg.eBOSC.trialPad          = 750;                                          % complete padding (WL + shoulder)
cfg.eBOSC.BGpad             = 750;                                          % padding of segments for BG (only avoiding edge artifacts)
cfg.eBOSC.fres              = 1.2;                                          % cf. Linkenkaer-Hansen, K., et al. (2001). "Long-Range Temporal Correlations and Scaling Behavior in Human Brain Oscillations." The Journal of Neuroscience 21(4): 1370-1377.
cfg.eBOSC.fstp              = 1;
cfg.eBOSC.freqRemoval       = 'JQK';
cfg.eBOSC.BiasCorrection    = 'yes';
cfg.eBOSC.LowFreqExcludeBG  = 8;
cfg.eBOSC.HighFreqExcludeBG = 15;

%%  initialize output matrices

load([pn.backgroundData 'background.mat'],'bckgrnd_filt')

amountAmps          = numel(cfg.simParams.amplitude);
amountCycles        = numel(cfg.simParams.cycles);
amountTimePoints    = numel(cfg.simParams.time);
amountFreqs         = numel(cfg.eBOSC.F);
amountRepetitions   = cfg.simParams.repetitions;

SignalDetection.Hits    = NaN(amountAmps, amountCycles, amountRepetitions,1);
SignalDetection.Misses  = NaN(amountAmps, amountCycles, amountRepetitions,1);
SignalDetection.CAs     = NaN(amountAmps, amountCycles, amountRepetitions,1);
SignalDetection.FAs     = NaN(amountAmps, amountCycles, amountRepetitions,1);

abundance_spec_ep   = NaN(amountAmps,amountCycles,amountRepetitions,amountFreqs);
abundance_ep        = NaN(amountAmps,amountCycles,amountRepetitions,1);
abundance_spec      = NaN(amountAmps,amountCycles,amountRepetitions, amountFreqs);

%%  loop conditions

for a = 1:length(cfg.simParams.amplitude)
    disp(num2str(a))
    for c = 1:length(cfg.simParams.cycles)
        disp(num2str(c))
        
        %% create data segments according to power and duration
        for k = 1:amountRepetitions
            % generate alpha in the middle of the segment
            rhythmCycles = cfg.simParams.cycles(c);
            rhythmFreq = 10;
            rhythmTime(c) = round((rhythmCycles/rhythmFreq),3);
            % simulate rhythms as symmetrical around the center
            timeNew = round(rhythmTime(c)/0.004,0);
            if mod(timeNew,2) ~= 0
                timeNew = timeNew + 1;
                rhythmTime(c) = timeNew.*0.004;
            else rhythmTime(c) = timeNew.*0.004;
            end; clear timeNew;
            rhythmTimeVector = [.004:.004:rhythmTime(c)];
            rhythmIdxVector = (amountTimePoints/2)-(numel(rhythmTimeVector)/2)+1:(amountTimePoints/2)+(numel(rhythmTimeVector)/2);
            % filter entire signal between 8 and 12 Hz (6th order butterworth) (not locally on alpha)
            [tmp_b,tmp_a] = butter(6, [8, 12]/(250/2), 'bandpass');
            tmp_bpsignal = filter(tmp_b,tmp_a,bckgrnd_filt(k,:));
            %VarBG = var(tmp_bpsignal(rhythmIdxVector));
            VarBG = var(tmp_bpsignal);
            % scale rhythm by noise background power
            targetPower = cfg.simParams.amplitude(a)*VarBG;
            amplitudeFromRMS = (sqrt(targetPower)*sqrt(2));
            simulatedRhythm = sin(rhythmTimeVector*2*pi*rhythmFreq)*amplitudeFromRMS;
            simulatedRhythm_complete = zeros(1,amountTimePoints);
            simulatedRhythm_complete(1,rhythmIdxVector) = simulatedRhythm;
            % effective segment = BG + rhythm
            signal = bckgrnd_filt(k,:) + simulatedRhythm_complete;
            %  wavelet transform to inspect 'simulated SNR'
            B(k,:,:) = BOSC_tf(signal,cfg.eBOSC.F,cfg.eBOSC.fsample,cfg.eBOSC.wavenumber);
            BGonly_FT(k,:,:) = BOSC_tf(bckgrnd_filt(k,:),cfg.eBOSC.F,cfg.eBOSC.fsample,cfg.eBOSC.wavenumber);
            freq = 19;
            FT_SNR_local(a,c,k) = mean(B(k,freq,rhythmIdxVector))./mean(BGonly_FT(k,freq,rhythmIdxVector));
            FT_SNR_total(a,c,k) = mean(B(k,19,rhythmIdxVector))./mean(BGonly_FT(k,19,:));
        end % k repetitions
        
        %%  run extended BOSC BG estimation
        
        % background is estimated across all simulated trials
        BG = [];
        for k = 1:amountRepetitions
            BG = cat(3,BG, B(k,:,:));
        end
        BG = squeeze(BG);
        
        % background power estimation - robust
        % find peak between 8-15 Hz, get wavelet extension in frequency domain, remove
        % points within this range from the estimation; compute robust
        % regression
        
        freqInd1 = find(cfg.eBOSC.F >= cfg.eBOSC.LowFreqExcludeBG, 1, 'first');
        freqInd2 = find(cfg.eBOSC.F <= cfg.eBOSC.HighFreqExcludeBG, 1, 'last');
        
        [~, indPos] = max(mean(BG(freqInd1:freqInd2,:),2));
        indPos = freqInd1+indPos;
        
        LowFreq = cfg.eBOSC.F(indPos)-(((2/cfg.eBOSC.wavenumber)*cfg.eBOSC.F(indPos))/2);
        UpFreq = cfg.eBOSC.F(indPos)+(((2/cfg.eBOSC.wavenumber)*cfg.eBOSC.F(indPos))/2);
        
        freqIndLow = find(cfg.eBOSC.F >= LowFreq, 1, 'first');
        freqIndHigh = find(cfg.eBOSC.F <= UpFreq, 1, 'last');
        
        [pv,~] = eBOSC_bgfit_robust(cfg.eBOSC.F([1:freqIndLow-1 freqIndHigh+1:end]),BG([1:freqIndLow-1 freqIndHigh+1:end], :));
        mp = 10.^(polyval(pv,log10(cfg.eBOSC.F)));
        
        % thresholds
        [pt,dt] = BOSC_thresholds(cfg.eBOSC.fsample,cfg.eBOSC.percentile,cfg.eBOSC.ncyc,cfg.eBOSC.F,mp);
        
        % keep overall background, 1/f fit, and power threshold
        BGinfo.all.bg_pow(a,c,k,:) = mean(BG(:,cfg.eBOSC.trialPad+1:end-cfg.eBOSC.trialPad),2);
        BGinfo.all.bg_log10_pow(a,c,k,:) = mean(log10(BG(:,cfg.eBOSC.trialPad+1:end-cfg.eBOSC.trialPad)),2);
        BGinfo.all.bg_amp(a,c,k,:) = mean(sqrt(BG(:,cfg.eBOSC.trialPad+1:end-cfg.eBOSC.trialPad)),2);
        BGinfo.all.pv(a,c,k,:) = pv;
        BGinfo.all.mp(a,c,k,:) = mp;
        BGinfo.all.pt(a,c,k,:)= pt;
        
        %% run rhythm detection across repetitions
        
        for k = 1:amountRepetitions
            
            display([num2str(a) '/8, ' num2str(c) '/8, ' num2str(k) '/500'])
            
            % oscillation detection
            detected = zeros(35,amountTimePoints-2*cfg.eBOSC.WLpadding);
            for f = 1:length(cfg.eBOSC.F)
                detected(f,:) = BOSC_detect(squeeze(B(k,f,cfg.eBOSC.WLpadding+1:end-cfg.eBOSC.WLpadding))',pt(f),dt(f),cfg.eBOSC.fsample);
            end; clear f
            
            % detected
            abundance_spec(a,c,k,:) = mean(detected(:,cfg.eBOSC.detectedPad+1:end-cfg.eBOSC.detectedPad),2);
            
            %%  abundance estimation
            
            cfg.eBOSC.npnts = size(detected,2);
            cfg.eBOSC.pt  = pt;
            
            if strcmp(variant, 'A')
                cfg.eBOSC.BiasVersion = '170630_MaxBias_all';
                cfg.eBOSC.BiasCorrection = 'yes';
                cfg.eBOSC.method = 'MaxBias';
                cfg.eBOSC.edgeOnly = 'no';
                cfg.eBOSC.effSignal = 'all';
            elseif strcmp(variant, 'B')
                cfg.eBOSC.BiasVersion = '170630_MaxBias_PT';
                cfg.eBOSC.BiasCorrection = 'yes';
                cfg.eBOSC.method = 'MaxBias';
                cfg.eBOSC.edgeOnly = 'no';
                cfg.eBOSC.effSignal = 'PT';
            elseif strcmp(variant, 'C')
                cfg.eBOSC.BiasVersion = '170630_HM_all';
                cfg.eBOSC.BiasCorrection = 'yes';
                cfg.eBOSC.method = 'HM';
                cfg.eBOSC.edgeOnly = 'no';
                cfg.eBOSC.effSignal = 'all';
            elseif strcmp(variant, 'D')
                cfg.eBOSC.BiasVersion = '170630_HM_PT';
                cfg.eBOSC.BiasCorrection = 'yes';
                cfg.eBOSC.method = 'HM';
                cfg.eBOSC.edgeOnly = 'no';
                cfg.eBOSC.effSignal = 'PT';
            elseif strcmp(variant, 'E')
                cfg.eBOSC.BiasCorrection = 'no';
            elseif strcmp(variant, 'G')
                cfg.eBOSC.BiasVersion = '171023_MaxBias_all';
                cfg.eBOSC.BiasCorrection = 'yes';
                cfg.eBOSC.method = 'MaxBias';
                cfg.eBOSC.edgeOnly = 'yes';
                cfg.eBOSC.effSignal = 'all';
            elseif strcmp(variant, 'H')
                cfg.eBOSC.BiasVersion = '171023_MaxBias_PT';
                cfg.eBOSC.BiasCorrection = 'yes';
                cfg.eBOSC.method = 'MaxBias';
                cfg.eBOSC.edgeOnly = 'yes';
                cfg.eBOSC.effSignal = 'PT';
            elseif strcmp(variant, 'I')
                cfg.eBOSC.BiasVersion = '171023_HM_all';
                cfg.eBOSC.BiasCorrection = 'yes';
                cfg.eBOSC.method = 'HM';
                cfg.eBOSC.edgeOnly = 'yes';
                cfg.eBOSC.effSignal = 'all';
            elseif strcmp(variant, 'J')
                cfg.eBOSC.BiasVersion = '171023_HM_PT';
                cfg.eBOSC.BiasCorrection = 'yes';
                cfg.eBOSC.method = 'HM';
                cfg.eBOSC.edgeOnly = 'yes';
                cfg.eBOSC.effSignal = 'PT';
            end
            
            episodes = [];
            [detected1,episodes] = eBOSC_createEpisodes(squeeze(B(k,:,cfg.eBOSC.WLpadding+1:end-cfg.eBOSC.WLpadding)),detected, cfg);
            
            abundance_spec_ep(a,c,k,:) = mean(detected1(:,cfg.eBOSC.detectedPad+1:end-cfg.eBOSC.detectedPad),2);
            
            segLength = size(detected1(:,cfg.eBOSC.detectedPad+1:end-cfg.eBOSC.detectedPad),2);
            
            % extra: plot detected vs detected1 power
            %figure; subplot(121); imagesc(detected.*squeeze(B(k,:,cfg.eBOSC.WLpadding+1:end-cfg.eBOSC.WLpadding))); subplot(122); imagesc(detected1.*squeeze(B(k,:,cfg.eBOSC.WLpadding+1:end-cfg.eBOSC.WLpadding))); 
            
            %% remove episodes and part of episodes that fall into 'shoulder'
            
            ind1 = cfg.eBOSC.detectedPad+1;
            ind2 = size(detected1,2) - cfg.eBOSC.detectedPad;
            cnt = 1;
            rmv = [];
            for j = 1:size(episodes,1)
                % final episodes
                ex = find(episodes{j,1}(:,2) < ind1 | episodes{j,1}(:,2) > ind2);
                % update episodes
                episodes{j,1}(ex,:) = [];
                episodes{j,1}(:,2) = episodes{j,1}(:,2) - ind1 + 1;
                episodes{j,2}(ex,:) = [];
                episodes{j,3}       = mean(episodes{j,2}(:,1));
                episodes{j,4}       = size(episodes{j,1},1) / cfg.eBOSC.fsample; clear ex
                if isempty(episodes{j,1})
                    rmv(cnt,1) = j;
                    cnt = cnt + 1;
                end
                % original episodes
                ex = find(episodes{j,5}(:,2) < ind1 | episodes{j,5}(:,2) > ind2);
                % update episodes
                episodes{j,5}(ex,:) = [];
                episodes{j,5}(:,2) = episodes{j,5}(:,2) - ind1 + 1;
                episodes{j,6}(ex,:) = [];
                episodes{j,7}       = mean(episodes{j,6}(:,1));
                episodes{j,8}       = size(episodes{j,5},1) / cfg.eBOSC.fsample; clear ex
            end; clear j cnt
            episodes(rmv,:) = []; clear rmv
            
            numberOfEffectivePoints = (amountTimePoints-2*cfg.eBOSC.trialPad);
            
            %% Signal Detection indices
            
            if ~isempty(episodes) % if episodes
                alpha_eps = find([episodes{:,3}]>= 8 & [episodes{:,3}]<= 12)';
                if ~isempty(alpha_eps) % if rhythmic episode
                    % create unique detected vector
                    alpha_locs = cat(1,episodes{alpha_eps,1});
                    alpha_locs = unique(alpha_locs(:,2))'; % detected rhythm SPs
                    % extract abundance
                    abundance_ep(a,c,k,1) = numel(alpha_locs)/numberOfEffectivePoints;
                    % check for overlap with simulated rhythm
                    if a > 1
                        % identify data points in final episode space that were
                        % simulated as rhythmic. Note that this gives an
                        % accurate amount of total signals as well.
                        alpha_sim_locs = zeros(1,amountTimePoints);
                        alpha_sim_locs(rhythmIdxVector) = 1;
                        alpha_sim_locs = find(alpha_sim_locs(cfg.eBOSC.trialPad+1:end-cfg.eBOSC.trialPad) == 1);
                    else
                        alpha_sim_locs = [];
                    end
                    
                    % encode amount of rhythm/non-rhythm
                    Amount.Alpha(c) = numel(alpha_sim_locs);
                    Amount.NoAlpha(c) = numberOfEffectivePoints-numel(alpha_sim_locs);
                    
                    Hits = intersect(alpha_locs, alpha_sim_locs); % detected & simulated
                    SignalDetection.Hits(a,c,k,1) = numel(Hits);
                    SignalDetection.HitRate(a,c,k,1) = numel(Hits)/numel(alpha_sim_locs);
                    FAs = setdiff(alpha_locs, alpha_sim_locs); % detected, not simulated
                    SignalDetection.FAs(a,c,k,1) = numel(FAs);
                    SignalDetection.FARate(a,c,k,1) = numel(FAs)/(numberOfEffectivePoints-numel(alpha_sim_locs));
                    Misses = setdiff(alpha_sim_locs, alpha_locs); % simulated, not detected
                    SignalDetection.Misses(a,c,k,1) = numel(Misses);
                    SignalDetection.MissRate(a,c,k,1) = numel(Misses)/numel(alpha_sim_locs);
                    CRs = setdiff(1:numberOfEffectivePoints, unique([alpha_sim_locs, alpha_locs])); % not simulated, not detected
                    SignalDetection.CRs(a,c,k,1) = numel(CRs);
                    SignalDetection.CRRate(a,c,k,1) = numel(CRs)/(numberOfEffectivePoints-numel(alpha_sim_locs));
                else
                    SignalDetection.FAs(a,c,k,1) = 0;
                    SignalDetection.FARate(a,c,k,1) = 0;
                    if a == 1
                        SignalDetection.Misses(a,c,k,1) = NaN;
                        SignalDetection.MissRate(a,c,k,1) = NaN;
                        SignalDetection.Hits(a,c,k,1) = NaN;
                        SignalDetection.HitRate(a,c,k,1) = NaN;
                        SignalDetection.CRs(a,c,k,1) = numberOfEffectivePoints;
                        SignalDetection.CRRate(a,c,k,1) = 1;
                    else
                        SignalDetection.Misses(a,c,k,1) = numel(alpha_sim_locs);
                        SignalDetection.MissRate(a,c,k,1) = 1;
                        SignalDetection.Hits(a,c,k,1) = 0;
                        SignalDetection.HitRate(a,c,k,1) = 0;
                        SignalDetection.CRs(a,c,k,1) = numberOfEffectivePoints-numel(alpha_sim_locs);
                        SignalDetection.CRRate(a,c,k,1) = 1;
                    end
                end
            else
                SignalDetection.FAs(a,c,k,1) = 0;
                SignalDetection.FARate(a,c,k,1) = 0;
                if a == 1
                    SignalDetection.Misses(a,c,k,1) = NaN;
                    SignalDetection.MissRate(a,c,k,1) = NaN;
                    SignalDetection.Hits(a,c,k,1) = NaN;
                    SignalDetection.HitRate(a,c,k,1) = NaN;
                    SignalDetection.CRs(a,c,k,1) = numberOfEffectivePoints;
                    SignalDetection.CRRate(a,c,k,1) = 1;
                else
                    SignalDetection.Misses(a,c,k,1) = numel(alpha_sim_locs);
                    SignalDetection.MissRate(a,c,k,1) = 1;
                    SignalDetection.Hits(a,c,k,1) = 0;
                    SignalDetection.HitRate(a,c,k,1) = 0;
                    SignalDetection.CRs(a,c,k,1) = numberOfEffectivePoints-numel(alpha_sim_locs);
                    SignalDetection.CRRate(a,c,k,1) = 1;
                end
            end
            
            %% extract rhythmic characteristics (e.g. power, frequency, etc.)
            
            if ~isempty(episodes)
                relevantEpisodeIDX = find([episodes{:,3}] >= 8 & [episodes{:,3}] <= 12);
                if isempty(relevantEpisodeIDX) || max(relevantEpisodeIDX) == 0 % if no episodes were found in the range
                    SignalDetection.Freq(a,c,k,1) = NaN;
                    SignalDetection.Amp(a,c,k,1) = NaN;
                    SignalDetection.oAmp(a,c,k,1) = NaN;
                    SignalDetection.fitBG(a,c,k,1) = NaN;
                    SignalDetection.PT(a,c,k,1) = NaN;
                    SignalDetection.BGdiff(a,c,k,1) = NaN;
                    SignalDetection.BGrel(a,c,k,1) = NaN;
                    SignalDetection.abn(a,c,k,1) = 0;
                    SignalDetection.abn_THG(a,c,k,1) = 0;
                    SignalDetection.NaNmat(a,c,k,1) = 1;
                else
                    %% frequency & amplitude
                    matchingEpisodes = cell2mat(episodes(relevantEpisodeIDX,2));
                    freqs = matchingEpisodes(:,1);
                    amps = sqrt(matchingEpisodes(:,2)); % square root of power values
                    %% row & column indices
                    Rows_Columns = cell2mat(episodes(relevantEpisodeIDX,1));
                    rows = Rows_Columns(:,1);
                    columns = Rows_Columns(:,2);
                    % unique column indices
                    uniqueSPs = unique(columns);
                    %% power spectrum of fitted BG
                    fitBG = sqrt(mp);
                    selectBG = fitBG(rows)';
                    % part of fitted BG that corresponds to frequencies of detected points
                    SignalDetection.fitBG(a,c,k,1) = nanmean(selectBG);
                    %% power threshold
                    SignalDetection.PT(a,c,k,1) = nanmean(pt(rows));
                    % absolute amplitude within episode
                    SignalDetection.Freq(a,c,k,1) = nanmean(freqs,1);
                    SignalDetection.Amp(a,c,k,1) = nanmean(amps,1);
                    % difference of detected amplitude from fitted BG
                    SignalDetection.BGdiff(a,c,k,1)  = nanmean(amps-selectBG,1);
                    SignalDetection.BGrel(a,c,k,1)   = nanmean((amps-selectBG)./selectBG,1);
                    SignalDetection.abn(a,c,k,1)     = numel(uniqueSPs)/segLength;
                    SignalDetection.abn_THG(a,c,k,1) = numel(Rows_Columns(:,2))/segLength;
                    SignalDetection.NaNmat(a,c,k,1)  = 0;
                end
            else
                SignalDetection.Freq(a,c,k,1) = NaN;
                SignalDetection.Amp(a,c,k,1) = NaN;
                SignalDetection.oAmp(a,c,k,1) = NaN;
                SignalDetection.fitBG(a,c,k,1) = NaN;
                SignalDetection.PT(a,c,k,1) = NaN;
                SignalDetection.BGdiff(a,c,k,1) = NaN;
                SignalDetection.BGrel(a,c,k,1) = NaN;
                SignalDetection.abn(a,c,k,1) = 0;
                SignalDetection.abn_THG(a,c,k,1) = 0;
                SignalDetection.NaNmat(a,c,k,1) = 1;
            end
            
            % overall amplitude of spectrum
            SpectrumOfInterest = squeeze(mean(B(k,17:21,cfg.eBOSC.trialPad+1:end-cfg.eBOSC.trialPad),3));
            mx = THG_find_local_max_1D_20150117(SpectrumOfInterest,1); % index of IAF peak
            if isempty(mx) & isnan(SpectrumOfInterest)
                SignalDetection.oAmp(a,c,k,1) = NaN;
                SignalDetection.oAmp_BGdiff(a,c,k,1) = NaN;
                SignalDetection.oAmp_BGrel(a,c,k,1) = NaN;
                SignalDetection.oBGestimate(a,c,k,1) = NaN;
                SignalDetection.oSNR(a,c,k,1) = NaN;
                SignalDetection.oFreq(a,c,k,1) = NaN;
            else
                if length(mx) == 1
                    mx = mx;
                elseif length(mx) > 1 % multiple peaks
                    mx = mx(SpectrumOfInterest(mx)==max(SpectrumOfInterest(mx))); % largest IAF peak
                elseif length(mx) < 1 % no peak
                    mx = find(SpectrumOfInterest==max(SpectrumOfInterest)); % largest IAF peak
                end
                % absolute amplitude
                SignalDetection.oAmp(a,c,k,1) = sqrt(SpectrumOfInterest(mx));
                % difference measure (from BG)
                fitBG = sqrt(mp);
                bgamp = fitBG(17+mx-1); % mean fitted BG spectrum at current channel & IAF
                SignalDetection.oAmp_BGdiff(a,c,k,1) = sqrt(SpectrumOfInterest(mx))-bgamp;
                SignalDetection.oAmp_BGrel(a,c,k,1) = (sqrt(SpectrumOfInterest(mx))-bgamp)./bgamp;
                SignalDetection.oBGestimate(a,c,k,1) = bgamp;
                SignalDetection.oSNR(a,c,k,1) = sqrt(SpectrumOfInterest(mx))./bgamp;
                SignalDetection.oFreq(a,c,k,1) = cfg.eBOSC.F(17+mx-1);
            end
            clear SpectrumOfInterest;
            
        end; clear k; % iterations
        
    end; clear c; % cycles
    
end; clear a; % amplitudes

save([pn.out, 'REDSimulation_171023_v10',variant,'.mat'], '-v7.3');

end