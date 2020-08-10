% Example script for eBOSC using simulations

% This script simulates alpha of varying cycles and amplitudes superimposed
% on pink noise and runs eBOSC on these simulated signals. Afterwards,
% signal detection measures are calculated.

%% set up paths
% automatically get script location from editor handle
tmp = matlab.desktop.editor.getActive;
pn.root = [fileparts(tmp.Filename), '/']; clear tmp;
addpath([pn.root, 'internal']) % add eBOSC functions
addpath([pn.root, 'external']) % add f_alpha_gaussian function
addpath([pn.root, 'external/BOSC']) % add BOSC functions

%% simulation parameters

cfg.simParams.amplitude     = [0 2 4 6 8 12 16 24];                         % simulated signal SNR
cfg.simParams.cycles        = [2 4 8 16 32 64 128 200];                     % simulated signal durations [in cycles]
cfg.simParams.segmentDur    = 20;                                           % duration of total segement [in seconds]
cfg.simParams.fsample       = 250;                                          % sampling rate of simulated signal [in Hz]
cfg.simParams.trials        = 100;                                          % amount of to-be-simulated trials
cfg.simParams.rhythmFreq    = 10;

%% create background data

% Note that the function 'f_alpha_gaussian' from the CNOISE toolbox 
% has to be added to the path. The toolbox is available at
% https://people.sc.fsu.edu/~jburkardt/m_src/cnoise/cnoise.html
% For reproducibility, we use a fixed seed.

% initialize background matrix
bckgrnd_filt = zeros(1000,5000);
% seed
randn('seed',20160118);
% loop 1000 repetitions
for k = 1:cfg.simParams.trials
    % generate 1/f background
    bckgrnd_filt(k,:) = f_alpha_gaussian(5000,1,1);
    % bandpass filter signal (consecutive low + high-pass filter)
    [B,A]  = butter(4,70/(250/2),'low'); 
    bckgrnd_filt(k,:) = filtfilt(B,A,bckgrnd_filt(k,:)); clear A B
    [B,A]  = butter(4,.5/(250/2),'high'); 
    bckgrnd_filt(k,:) = filtfilt(B,A,bckgrnd_filt(k,:)); clear A B
end; clear k

%% create simulated data structure with specified power and duration

count = 1;
for a = 1:length(cfg.simParams.amplitude)
    for c = 1:length(cfg.simParams.cycles)
        for k = 1:cfg.simParams.trials
            data.time{k} = [1/cfg.simParams.fsample:1/cfg.simParams.fsample:cfg.simParams.segmentDur];
            % generate alpha in the middle of the segment
            rhythmTime(c) = round((cfg.simParams.cycles(c)/cfg.simParams.rhythmFreq),3);
            % simulate rhythms as symmetrical around the center
            timeNew = round(rhythmTime(c)/0.004,0);
            if mod(timeNew,2) ~= 0
                timeNew = timeNew + 1;
                rhythmTime(c) = timeNew.*0.004;
            else rhythmTime(c) = timeNew.*0.004;
            end; clear timeNew;
            rhythmTimeVector = [.004:.004:rhythmTime(c)];
            rhythmIdxVector = (numel(data.time{k})/2)-(numel(rhythmTimeVector)/2)+1:(numel(data.time{k})/2)+(numel(rhythmTimeVector)/2);
            % filter entire signal between 8 and 12 Hz (6th order butterworth) (not locally on alpha)
            [tmp_b,tmp_a] = butter(6, [8, 12]/(250/2), 'bandpass');
            tmp_bpsignal = filter(tmp_b,tmp_a,squeeze(bckgrnd_filt(k,:))); clear tmp_a tmp_b
            % scale rhythm by noise background power
            amplitudeFromRMS = (sqrt(cfg.simParams.amplitude(a)*var(tmp_bpsignal))*sqrt(2)); clear tmp_bpsignal;
            simulatedRhythm = sin(rhythmTimeVector*2*pi*cfg.simParams.rhythmFreq)*amplitudeFromRMS;
            simulatedRhythm_complete = zeros(1,numel(data.time{k}));
            simulatedRhythm_complete(1,rhythmIdxVector) = simulatedRhythm;
            % effective segment = BG + rhythm
            data.trial{k}(count,:) = bckgrnd_filt(k,:) + simulatedRhythm_complete;
            % create an artificial channel with current amplitude x
            % duration combination
            data.label{count} = ['a_', num2str(a), '_c_', num2str(c)];
            % encode when the rhythm was simulated
            data.rhythm(count,:) = zeros(size(simulatedRhythm_complete));
            data.rhythm(count,rhythmIdxVector) = 1;
        end; clear simulatedRhythm_complete simulatedRhythm amplitudeFromRMS
        count = count + 1;
    end
end; clear count rhythmTimeVector rhythmIdxVector rhythmTime a c k bckgrnd_filt

%% eBOSC parameters

% general setup
cfg.eBOSC.F             = 2.^[1:.125:5.25];
cfg.eBOSC.wavenumber	= 6;
cfg.eBOSC.fsample       = cfg.simParams.fsample;

% padding
cfg.eBOSC.pad.tfr_s = 2;           % padding following wavelet transform to avoid edge artifacts in seconds (bi-lateral)
cfg.eBOSC.pad.detection_s = 1;     % padding following rhythm detection in seconds (bi-lateral); 'shoulder' for BOSC eBOSC.detected matrix to account for duration threshold
cfg.eBOSC.pad.background_s = 2;    % padding of segments for BG (only avoiding edge artifacts)

% threshold settings
cfg.eBOSC.threshold.excludePeak = [8,15];                                   % lower and upper bound of frequencies to be excluded during background fit (Hz) (previously: LowFreqExcludeBG HighFreqExcludeBG)
cfg.eBOSC.threshold.duration	= repmat(3, 1, numel(cfg.eBOSC.F));         % vector of duration thresholds at each frequency (previously: ncyc)
cfg.eBOSC.threshold.percentile  = .95;                                      % percentile of background fit for power threshold

% episode post-processing
cfg.eBOSC.postproc.use      = 'no';         % Post-processing of rhythmic eBOSC.episodes, i.e., wavelet 'deconvolution' (default = 'no')
cfg.eBOSC.postproc.method   = 'MaxBias';	% Deconvolution method (default = 'MaxBias', FWHM: 'FWHM')
cfg.eBOSC.postproc.edgeOnly = 'yes';        % Deconvolution only at on- and offsets of eBOSC.episodes? (default = 'yes')
cfg.eBOSC.postproc.effSignal= 'PT';         % Power deconvolution on whole signal or signal above power threshold? (default = 'PT')

% general processing settings
cfg.eBOSC.channel = []; % select channels (default: all)
cfg.eBOSC.trial = []; % select trials (default: all)
cfg.eBOSC.trial_background = []; % select trials for background (default: all)

%% run eBOSC

[eBOSC, cfg] = eBOSC_wrapper(cfg, data);

%% Sanity-check plots

indChan = 62; indTrial = 1; % Here we select the first trial and 62nd (artificial) channel we encoded (see cfg.eBOSC.channel).

disp(['Results are for trial ', num2str(cfg.eBOSC.trial(indTrial)), ' at channel ', data.label{cfg.eBOSC.channel(indChan)}])

% get original time series for plotting
origData = data.trial{indTrial}(cfg.eBOSC.channel(indChan), cfg.eBOSC.pad.total_sample+1:end-cfg.eBOSC.pad.total_sample);

% Supplementary Figure: plot estimated background + power threshold
figure; hold on;
plot(log10(cfg.eBOSC.F), log10(eBOSC.static.mp(indChan,:)), 'k--','LineWidth', 1.5); 
plot(log10(cfg.eBOSC.F), log10(eBOSC.static.pt(indChan,:)), 'k-', 'LineWidth', 1.5)
plot(log10(cfg.eBOSC.F), eBOSC.static.bg_log10_pow(indChan,:), 'r-', 'LineWidth', 2)
xlabel('Frequency (log10 Hz)'); ylabel('Power (log 10 a.u.)');
legend({'Aperiodic fit', 'Statistical power threshold', 'Avg. spectrum'}, ...
    'orientation', 'vertical', 'location', 'SouthWest'); legend('boxoff');
set(findall(gcf,'-property','FontSize'),'FontSize',20)
xlim([.3, 1.75])

% Supplementary Plot: eBOSC's detected matrix: plot rhythmic alpha episodes
h = figure('units','normalized','position',[.1 .1 .6 .3]);
hold on; 
plot(squeeze(origData), 'k');
tmpDetected = single(squeeze(nanmean(eBOSC.detected(indChan, indTrial,cfg.eBOSC.F > 8 & cfg.eBOSC.F < 15,:),3))>0); tmpDetected(tmpDetected==0) = NaN;
plot(squeeze(origData).*tmpDetected', 'r');
xlabel('Time (s)'); ylabel('Power [µV]');
legend({'Original signal'; 'Rhythmic signal'}, ...
    'orientation', 'horizontal', 'location', 'north'); legend('boxoff')
set(findall(gcf,'-property','FontSize'),'FontSize',26)

% Supplementary Plots: different episode statistics
figure; 
subplot(3,2,1); histogram(eBOSC.episodes.SNRMean); title('SNR distribution')
subplot(3,2,2); histogram(log10(eBOSC.episodes.SNRMean)); title('SNR distribution(log10)')
subplot(3,2,3); histogram(eBOSC.episodes.DurationC); title('Duration distribution')
subplot(3,2,4); histogram(log10(eBOSC.episodes.DurationC)); title('Duration distribution(log10)')
subplot(3,2,5); histogram(eBOSC.episodes.FrequencyMean); title('Frequency distribution')
subplot(3,2,6); hold on; plot(squeeze(eBOSC.pepisode(indChan, indTrial,:))); plot(squeeze(eBOSC.abundance_ep(indChan, indTrial,:))); title('Pepisode, abundance')

% Supplementary Plot: plot rhythmic episodes with indicated onsets
idx_alpha = find(eBOSC.episodes.Trial == indTrial & eBOSC.episodes.Channel == cfg.eBOSC.channel(indChan) &...
    eBOSC.episodes.FrequencyMean > 8 & eBOSC.episodes.FrequencyMean <15);
% filter for alpha by mean frequency (!) of episode
idx_onset = []; idx_onsetTime = [];
alphaDetected = NaN(1,numel(origData));
for indEp = 1:numel(idx_alpha)
    % These are two alternative ways to extract the onset timepoint from the table
    idx_onsetTime(indEp) = find(cfg.tmp.finalTime>= eBOSC.episodes.Onset(idx_alpha(indEp)), 1, 'first');
    idx_onset(indEp) = eBOSC.episodes.ColID{idx_alpha(indEp)}(1);
    % Mark all periods with episodes falling into the alpha range
    alphaDetected(eBOSC.episodes.ColID{idx_alpha(indEp)}(1):eBOSC.episodes.ColID{idx_alpha(indEp)}(end)) = 1;
end; clear idx_alpha;
h = figure('units','normalized','position',[.1 .1 .7 .3]); hold on; 
scatter(idx_onset, repmat(100,1,numel(idx_onset)), 75, [.5 .5 .5], 'filled')
OnsetLine = squeeze(origData);
OnsetLine(idx_onset) = 100; clear idx_onset idx_onsetTime;
plot(OnsetLine, 'Color', [.5 .5 .5]); clear OnsetLine;
[orig]=plot(squeeze(origData), 'k');
[rhythm]=plot(squeeze(origData).*alphaDetected, 'r');
xlabel('Time (s)'); ylabel('Power [µV]');
legend([orig, rhythm], {'Original signal'; 'Rhythmic signal'}, ...
    'orientation', 'horizontal', 'location', 'south'); legend('boxoff')
set(findall(gcf,'-property','FontSize'),'FontSize',26)

%% Signal Detection indices

N_amp = numel(cfg.simParams.amplitude);
N_cyc = numel(cfg.simParams.cycles);
N_trial = cfg.simParams.trials;

% initialize with NaN
SignalDetection.MissRate = NaN(N_amp, N_cyc, N_trial, 1);
SignalDetection.HitRate = NaN(N_amp, N_cyc, N_trial, 1);
SignalDetection.CRRate = NaN(N_amp, N_cyc, N_trial, 1);
SignalDetection.FARate = NaN(N_amp, N_cyc, N_trial, 1);

for a = 1:N_amp
    for c = 1:N_cyc
        for k = 1:N_trial
            % get relevant channel
            channelIdx = find(ismember(data.label, ['a_', num2str(a), '_c_', num2str(c)]));
            if isempty(channelIdx)
                error('Channel not found');
                continue;
            end
            % get relevant detected episodes
            alpha_eps = find(eBOSC.episodes.FrequencyMean>= 8 & eBOSC.episodes.FrequencyMean<= 12 & ...
                eBOSC.episodes.Trial== k & eBOSC.episodes.Channel== channelIdx)';
            if ~isempty(alpha_eps) % if a rhythmic episode with the criteria was found, encode results
                % create vector of detected timoints for episodes with criteria
                alpha_locs = cat(2,eBOSC.episodes.ColID{alpha_eps,1});
                alpha_locs = cellfun(@(x,y)x:y,num2cell(alpha_locs(1,:)),num2cell(alpha_locs(2,:)),'UniformOutput',false);
                alpha_locs = cat(2,alpha_locs{:});
                alpha_locs = unique(alpha_locs)';
                % get the number of eligible time points after accounting for padding
                numberOfEffectivePoints = size(data.time{1},2) - 2*cfg.eBOSC.pad.total_sample;
                % check for overlap with simulated rhythm
                if cfg.simParams.amplitude(a) > 0
                    % identify data points in final episode space that were
                    % simulated as rhythmic. Note that this gives an
                    % accurate amount of total signals as well.
                    alpha_sim_locs = find(data.rhythm(channelIdx,cfg.eBOSC.pad.total_sample+1:end-cfg.eBOSC.pad.total_sample));
                else
                    % This condition captures the case where the rhythmic
                    % amplitude was 0.
                    alpha_sim_locs = [];
                end
                Hits = intersect(alpha_locs, alpha_sim_locs); % detected & simulated
                FAs = setdiff(alpha_locs, alpha_sim_locs); % detected, not simulated
                Misses = setdiff(alpha_sim_locs, alpha_locs); % simulated, not detected
                CRs = setdiff(1:numberOfEffectivePoints, unique([alpha_sim_locs, alpha_locs'])); % not simulated, not detected
                
                SignalDetection.HitRate(a,c,k,1) = numel(Hits)/numel(alpha_sim_locs);
                SignalDetection.FARate(a,c,k,1) = numel(FAs)/(numberOfEffectivePoints-numel(alpha_sim_locs));
                SignalDetection.MissRate(a,c,k,1) = numel(Misses)/numel(alpha_sim_locs);
                SignalDetection.CRRate(a,c,k,1) = numel(CRs)/(numberOfEffectivePoints-numel(alpha_sim_locs));
            else
                SignalDetection.FARate(a,c,k,1) = 0;
                if cfg.simParams.amplitude(a) == 0
                    SignalDetection.MissRate(a,c,k,1) = NaN;
                    SignalDetection.HitRate(a,c,k,1) = NaN;
                    SignalDetection.CRRate(a,c,k,1) = 1;
                else
                    SignalDetection.MissRate(a,c,k,1) = 1;
                    SignalDetection.HitRate(a,c,k,1) = 0;
                    SignalDetection.CRRate(a,c,k,1) = 1;
                end
            end
        end; clear k; % iterations
    end
end

% plot signal detection overview
figure; hold on;
scatter(squeeze(SignalDetection.FARate(:)),squeeze(SignalDetection.HitRate(:)),10, 'k', 'filled')
xlabel('False positive rate') 
ylabel('True positive rate')
