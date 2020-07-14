% Example script for eBOSC analysis

%% initial clean-up

clear all; clc; restoredefaultpath;

%% add eBOSC toolbox with subdirectories

% automatically get script location from editor handle
tmp = matlab.desktop.editor.getActive;
pn.eBOSC = [fileparts(tmp.Filename), '/']; clear tmp;

cd(pn.eBOSC)
addpath(pn.eBOSC);
addpath([pn.eBOSC, 'dev/']);
addpath([pn.eBOSC, 'external/BOSC/']);
% addpath([pn.eBOSC, 'external/fieldtrip/']); ft_defaults;
% addpath([pn.eBOSC, 'external/NoiseTools/']);

%% eBOSC parameters

% general setup
cfg.eBOSC.F             = 2.^[1:.125:6];    % frequency sampling (~Whitten et al., 2011), but higher frequency resolution
cfg.eBOSC.wavenumber	= 6;                % wavelet family parameter (time-frequency tradeoff) [recommended: ~6]
cfg.eBOSC.fsample       = 500;              % current sampling frequency of EEG data

% padding
cfg.eBOSC.pad.tfr_s = 1;                                                                    % padding following wavelet transform to avoid edge artifacts in seconds (bi-lateral)
cfg.eBOSC.pad.tfr_sample = cfg.eBOSC.pad.tfr_s.*cfg.eBOSC.fsample;                          % automatic sample rate calculation
cfg.eBOSC.pad.detection_s = .5;                                                             % padding following rhythm detection in seconds (bi-lateral); 'shoulder' for BOSC eBOSC.detected matrix to account for duration threshold
cfg.eBOSC.pad.detection_sample = cfg.eBOSC.pad.detection_s.*cfg.eBOSC.fsample;              % automatic sample rate calculation
cfg.eBOSC.pad.total_s = cfg.eBOSC.pad.tfr_s + cfg.eBOSC.pad.detection_s;                    % complete padding (WL + shoulder)
cfg.eBOSC.pad.total_sample = cfg.eBOSC.pad.tfr_sample + cfg.eBOSC.pad.detection_sample;
cfg.eBOSC.pad.background_s = cfg.eBOSC.pad.tfr_s;                                           % padding of segments for BG (only avoiding edge artifacts)
cfg.eBOSC.pad.background_sample = cfg.eBOSC.pad.tfr_sample;

% threshold settings
cfg.eBOSC.threshold.excludePeak = [8,15];                                   % lower and upper bound of frequencies to be excluded during background fit (Hz) (previously: LowFreqExcludeBG HighFreqExcludeBG)
cfg.eBOSC.threshold.duration	= repmat(3, 1, numel(cfg.eBOSC.F));         % vector of duration thresholds at each frequency (previously: ncyc)
cfg.eBOSC.threshold.percentile  = .95;                                      % percentile of background fit for power threshold

% episode post-processing
cfg.eBOSC.postproc.use      = 'no';         % Post-processing of rhythmic eBOSC.episodes, i.e., wavelet 'deconvolution' (default = 'no')
cfg.eBOSC.postproc.method   = 'MaxBias';	% Deconvolution method (default = 'MaxBias', FWHM: 'FWHM')
cfg.eBOSC.postproc.edgeOnly = 'yes';        % Deconvolution only at on- and offsets of eBOSC.episodes? (default = 'yes')
cfg.eBOSC.postproc.effSignal= 'PT';         % Amplitude deconvolution on whole signal or signal above power threshold? (default = 'PT')

%% load data

load([pn.eBOSC,  'util/1160_rest_EEG_Rlm_Fhl_rdSeg_Art_EC.mat'], 'data')

%% concatenate trials for resting state here

data.trial{1} = cat(2,data.trial{:}); data.trial(2:end) = [];
data.time{1} = cat(2,data.time{:}); data.time(2:end) = [];
data = rmfield(data, 'sampleinfo');

%% initialize eBOSC output structure

eBOSC = [];

%% ---- select a channel here

e = 20;
display(['channel #' num2str(e)])

cfg.tmp.inputTime = data.time{1,1};
cfg.tmp.detectedTime = cfg.tmp.inputTime(cfg.eBOSC.pad.tfr_sample+1:end-cfg.eBOSC.pad.tfr_sample);
cfg.tmp.finalTime = cfg.tmp.inputTime(cfg.eBOSC.pad.total_sample+1:end-cfg.eBOSC.pad.total_sample);
cfg.tmp.channel = e; % encode current channel for later

%% Step 1: time-frequency wavelet decomposition for whole signal to prepare background fit

eBOSC.Ntrial = length(data.trial);

TFR = [];
for indTrial = 1:eBOSC.Ntrial
    % get data
    tmp_dat = data.trial{indTrial}(e,:);
    % wavelet transform (NOTE: no check to avoid spectral leakage);
    % apply correction factor
    TFR.trial{indTrial} = BOSC_tf(tmp_dat,cfg.eBOSC.F,cfg.eBOSC.fsample,cfg.eBOSC.wavenumber);
    clear tmp_dat
end; clear indTrial

%% Step 2: robust background power fit (see 2020 NeuroImage paper)

[eBOSC, pt, dt] = eBOSC_getThresholds(cfg, TFR, eBOSC);

% Supplementary Figure: plot estimated background + power threshold
figure; hold on;
plot(log10(cfg.eBOSC.F), log10(eBOSC.static.mp(e,:)), 'k--','LineWidth', 1.5); 
plot(log10(cfg.eBOSC.F), log10(eBOSC.static.pt(e,:)), 'k-', 'LineWidth', 1.5)
plot(log10(cfg.eBOSC.F),log10(eBOSC.static.bg_pow(e,:)), 'r-', 'LineWidth', 2)
xlabel('Frequency (log10 Hz)'); ylabel('Power (log 10 a.u.)');
legend({'Aperiodic fit', 'Statistical power threshold', 'Avg. spectrum'}, ...
    'orientation', 'vertical', 'location', 'SouthWest'); legend('boxoff');
set(findall(gcf,'-property','FontSize'),'FontSize',20)
xlim([.3, 1.75])

%% application of thresholds to single trials

for indTrial = 1:eBOSC.Ntrial
        
    cfg.tmp.trial = indTrial; % encode current trial for later

    % get wavelet transform for single trial
    % WLpadding is removed to avoid edge artifacts during the
    % detection. Note that detectedPad still remains so that there
    % is no problems with too few sample points at the edges to
    % fulfill the numcycles criterion.
    TFR_ = TFR.trial{1,indTrial}(:,cfg.eBOSC.pad.tfr_sample+1:end-cfg.eBOSC.pad.tfr_sample);

    %% Step 3: detect rhythms and calculate Pepisode

    % The next section applies both the power and the optional duration
    % threshold to detect individual rhythmic segments in the continuous signals.

    detected = zeros(size(TFR_));
    for f = 1:length(cfg.eBOSC.F)
        detected(f,:) = BOSC_detect(TFR_(f,:),pt(f),dt(f),cfg.eBOSC.fsample);
    end; clear f

    % remove padding for detection (matrix with padding required for refinement)
    eBOSC.detected(indTrial,:,:) = detected(:,cfg.eBOSC.pad.detection_sample+1:end-cfg.eBOSC.pad.detection_sample);
    curDetected = squeeze(eBOSC.detected(indTrial,:,:));
    
    % encode pepisode of detected rhythms (optional)
    eBOSC.pepisode(indTrial,:) = mean(curDetected,2);
    
    % encode detected alpha signals (optional)
    alphaDetected = zeros(1,size(curDetected,2));
    alphaDetected(nanmean(curDetected(cfg.eBOSC.F > 8 & cfg.eBOSC.F < 15,:),1)>0) = 1;
    eBOSC.detectedAlpha(indTrial,:) = alphaDetected;

    % encode original signals (optional)
    origData = data.trial{indTrial}(e, cfg.eBOSC.pad.total_sample+1:end-cfg.eBOSC.pad.total_sample);
    eBOSC.origData(indTrial,:) = origData;

    % Supplementary Plot: plot only rhythmic episodes
    h = figure('units','normalized','position',[.1 .1 .6 .3]);
    hold on; 
    plot(eBOSC.origData(indTrial,:), 'k');
    tmpDetected = eBOSC.detectedAlpha(indTrial,:); tmpDetected(tmpDetected==0) = NaN;
    plot(eBOSC.origData(indTrial,:).*tmpDetected, 'r');
    xlim([7.2, 7.9]*10^4)
    xlabel('Time (s)'); ylabel('Amplitude [µV]');
    legend({'Original signal'; 'Rhythmic signal'}, ...
        'orientation', 'horizontal', 'location', 'north'); legend('boxoff')
    set(findall(gcf,'-property','FontSize'),'FontSize',26)
    
    %% Step 4 (optional): create table of separate rhythmic episodes

    eBOSC.episodes = [];
    [detected_ep,eBOSC.episodes] = eBOSC_episode_create(TFR_,eBOSC, cfg, detected);

    % remove padding for detection (already done for eBOSC.episodes)
    eBOSC.detected_ep(indTrial,:,:) = detected_ep(:,cfg.eBOSC.pad.detection_sample+1:end-cfg.eBOSC.pad.detection_sample);
    clear detected_ep;
    
    % encode abundance of eBOSC.episodes (optional)
    eBOSC.abundance_ep(indTrial,:) = mean(squeeze(eBOSC.detected_ep(indTrial,:,:)),2);
    
    % Supplementary Plot: original eBOSC.detected vs. sparse episode power
    figure; 
    subplot(121); imagesc(squeeze(eBOSC.detected).*TFR_(:,cfg.eBOSC.pad.detection_sample+1:end-cfg.eBOSC.pad.detection_sample));
    subplot(122); imagesc(squeeze(eBOSC.detected_ep).*TFR_(:,cfg.eBOSC.pad.detection_sample+1:end-cfg.eBOSC.pad.detection_sample));
    
    % Supplementary Plots: different statistics
    figure; 
    subplot(3,2,1); histogram(eBOSC.episodes.SNRMean); title('SNR distribution')
    subplot(3,2,2); histogram(log10(eBOSC.episodes.SNRMean)); title('SNR distribution(log10)')
    subplot(3,2,3); histogram(eBOSC.episodes.DurationC); title('Duration distribution')
    subplot(3,2,4); histogram(log10(eBOSC.episodes.DurationC)); title('Duration distribution(log10)')
    subplot(3,2,5); histogram(eBOSC.episodes.FrequencyMean); title('Frequency distribution')
    subplot(3,2,6); hold on; plot(squeeze(eBOSC.pepisode(indTrial,:))); plot(squeeze(eBOSC.abundance_ep(indTrial,:))); title('Pepisode, abundance')
    
    % encode detected alpha signals (optional)
    curDetected = squeeze(eBOSC.detected_ep(indTrial,:,:));
    alphaDetected = zeros(1,size(curDetected,2));
    alphaDetected(nanmean(curDetected(cfg.eBOSC.F > 8 & cfg.eBOSC.F < 15,:),1)>0) = 1;
    eBOSC.detectedAlpha_ep(indTrial,:) = alphaDetected;

    % Supplementary Plot: plot only rhythmic episodes
    figure; hold on; 
    plot(eBOSC.origData(indTrial,:), 'k');
    tmpDetected = eBOSC.detectedAlpha_ep(indTrial,:); tmpDetected(tmpDetected==0) = NaN;
    plot(eBOSC.origData(indTrial,:).*tmpDetected, 'r');
    xlim([7.2, 7.9]*10^4)
    
    % plot example of onsets for alpha signals
    
    % filter for alpha
    
    idx_alpha = find(eBOSC.episodes.FrequencyMean > 8 & eBOSC.episodes.FrequencyMean <15);
    
    exampleAlphaOnsetTime = []; exampleAlphaOnset = [];
    for indEp = 1:numel(idx_alpha)
        idx_onsetTime(indEp) = find(cfg.tmp.finalTime>= eBOSC.episodes.Onset(idx_alpha(indEp)), 1, 'first');
        idx_onset(indEp) = eBOSC.episodes.ColID{idx_alpha(indEp)}(1);
        % encode in matrix
        %exampleAlphaOnsetTime(indEp, :) = data.trial{indTrial}(e,idx_onsetTime(indEp)-500:idx_onsetTime(indEp)+500);
        %exampleAlphaOnset(indEp, :) = eBOSC.origData(indTrial,idx_onset:idx_onset+500);
    end

    figure; imagesc(exampleAlphaOnset)
    figure; plot(exampleAlphaOnset')

    % Supplementary Plot: plot only rhythmic episodes
    figure; hold on; 
    plot(eBOSC.origData(indTrial,:), 'k');
    tmpDetected = eBOSC.detectedAlpha_ep(indTrial,:); tmpDetected(tmpDetected==0) = NaN;
    plot(eBOSC.origData(indTrial,:).*tmpDetected, 'r');
    scatter(idx_onset, repmat(100,1,numel(idx_onset)), 'filled')
    OnsetLine = zeros(size(eBOSC.origData(indTrial,:)));
    OnsetLine(idx_onset) = 100;
    plot(OnsetLine, 'g')
    xlim([7.2, 7.9]*10^4)
    
    
end; clear indTrial;
