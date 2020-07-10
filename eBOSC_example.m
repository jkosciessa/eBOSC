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
addpath([pn.eBOSC, 'external/fieldtrip/']); ft_defaults;
addpath([pn.eBOSC, 'external/NoiseTools/']);

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

% episode creation
cfg.eBOSC.fstp = 1;

% episode post-processing
cfg.eBOSC.postproc.use      = 'yes';    % Post-processing of rhythmic eBOSC.episodes, i.e., wavelet 'deconvolution' (default = 'no')
cfg.eBOSC.postproc.method   = 'FWHM';	% Deconvolution method (default = 'MaxBias', FWHM: 'FWHM')
cfg.eBOSC.postproc.edgeOnly = 'yes';	% Deconvolution only at on- and offsets of eBOSC.episodes? (default = 'yes')
cfg.eBOSC.postproc.effSignal= 'PT';     % Amplitude deconvolution on whole signal or signal above power threshold? (default = 'PT')

%% load data

load([pn.eBOSC,  'util/1160_rest_EEG_Rlm_Fhl_rdSeg_Art_EC.mat'], 'data')

%% concatenate trials for resting state here

data.trial{1} = cat(2,data.trial{:}); data.trial(2:end) = [];
data.time{1} = cat(2,data.time{:}); data.time(2:end) = [];
data = rmfield(data, 'sampleinfo');

%% ---- select a channel here

e = 60;
display(['channel #' num2str(e)])

%% TO DO: include the following in 'dedicated' eBOSC function

% initialize eBOSC output structure

eBOSC = [];

% -----------------------------
% TO DO: include fieldtrip-style automatic detection of dimord
% -----------------------------

cfg.tmp.inputTime = data.time{1,1};
cfg.tmp.detectedTime = cfg.tmp.inputTime(cfg.eBOSC.pad.tfr_sample+1:end-cfg.eBOSC.pad.tfr_sample);
cfg.tmp.finalTime = cfg.tmp.inputTime(cfg.eBOSC.pad.total_sample+1:end-cfg.eBOSC.pad.total_sample);
cfg.tmp.channel = e; % encode current channel for later

%% TF analysis for whole signal to prepare background fit

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

%% eBOSC background: robust background power fit (2020 NeuroImage paper)

[eBOSC, pt, dt] = eBOSC_getThresholds(cfg, TFR, eBOSC);

% Supplementary Figure: plot estimated background + power threshold
figure; hold on;
plot(log10(cfg.eBOSC.F), log10(eBOSC.static.mp(e,:)), 'k--','LineWidth', 1.5); 
plot(log10(cfg.eBOSC.F), log10(eBOSC.static.pt(e,:)), 'k-', 'LineWidth', 1.5)
plot(log10(cfg.eBOSC.F),log10(eBOSC.static.bg_pow(e,:)), 'r-', 'LineWidth', 2)
xlabel('Frequency (log10 Hz)'); ylabel('Power (log 10 a.u.)');
legend({'Aperiodic fit', 'Statistical power threshold', 'Avg. spectrum'}, ...
    'orientation', 'vertical', 'location', 'SouthWest'); legend('boxoff');
set(findall(gcf,'-property','FontSize'),'FontSize',14)

%% use statically-defined threshold for single-trial detection

for indTrial = 1:eBOSC.Ntrial
        
    cfg.tmp.trial = indTrial; % encode current channel for later

    % get wavelet transform for single trial
    % WLpadding is removed to avoid edge artifacts during the
    % detection. Note that detectedPad still remains so that there
    % is no problems with too few sample points at the edges to
    % fulfill the numcycles criterion.
    TFR_ = TFR.trial{1,indTrial}(:,cfg.eBOSC.pad.tfr_sample+1:end-cfg.eBOSC.pad.tfr_sample);

%%  Detect rhythms + calculate Pepisode (~standard BOSC)

    % The next section applies both the power and the optional duration
    % threshold to detect individual rhythmic segments in the continuous signals.

    eBOSC.detected = zeros(size(TFR_));
    for f = 1:length(cfg.eBOSC.F)
        eBOSC.detected(f,:) = BOSC_detect(TFR_(f,:),pt(f),dt(f),cfg.eBOSC.fsample);
    end; clear f

    % encode pepisode of detected rhythms (optional)
    eBOSC.pepisode(indTrial,:) = mean(eBOSC.detected(:,cfg.eBOSC.pad.detection_sample+1:end-cfg.eBOSC.pad.detection_sample),2);
    
    % encode eBOSC.detected alpha signals (optional)
    alphaDetected = zeros(1,size(eBOSC.detected,2));
    alphaDetected(nanmean(eBOSC.detected(cfg.eBOSC.F > 8 & cfg.eBOSC.F < 12,:),1)>0) = 1;
    eBOSC.detectedAlpha(e,:) = alphaDetected(1,cfg.eBOSC.pad.detection_sample+1:end-cfg.eBOSC.pad.detection_sample);

    % encode original signals (optional)
    origData = data.trial{1}(e, cfg.eBOSC.pad.total_sample+1:end-cfg.eBOSC.pad.total_sample);
    eBOSC.origData(e,:) = origData;

    % Supplementary Plot: plot only rhythmic eBOSC.episodes
    figure; hold on; 
    plot(eBOSC.origData(e,:));
    tmpDetected = eBOSC.detectedAlpha; tmpDetected(tmpDetected==0) = NaN;
    plot(eBOSC.origData(e,:).*tmpDetected(e,:));
    xlim([7.2, 7.9]*10^4)

    %% create list of rhythmic eBOSC.episodes

    eBOSC.episodes = [];
    [detected_ep,eBOSC.episodes] = eBOSC_episode_create(TFR_,eBOSC, cfg);

    % encode abundance of eBOSC.episodes (optional)
    eBOSC.abundance_ep(indTrial,:) = mean(detected_ep(:,cfg.eBOSC.pad.detection_sample+1:end-cfg.eBOSC.pad.detection_sample),2);

    % Supplementary Plot: original eBOSC.detected vs. sparse episode power
    figure; 
    subplot(121); imagesc(eBOSC.detected.*TFR_);
    subplot(122); imagesc(eBOSC.detected_ep.*TFR_);
    
    % remove padding for detection (already done for eBOSC.episodes)
    eBOSC.detected_ep(indTrial,:) = detected_ep(:,cfg.eBOSC.pad.detection_sample+1:end-cfg.eBOSC.pad.detection_sample);
    clear detected_ep;

    % Supplementary Plots:

    figure; 
    subplot(3,2,1); histogram(eBOSC.episodes.SNRMean); title('SNR distribution')
    subplot(3,2,2); histogram(log10(eBOSC.episodes.SNRMean)); title('SNR distribution(log10)')
    subplot(3,2,3); histogram(eBOSC.episodes.DurationC); title('Duration distribution')
    subplot(3,2,4); histogram(log10(eBOSC.episodes.DurationC)); title('Duration distribution(log10)')
    subplot(3,2,5); histogram(eBOSC.episodes.FrequencyMean); title('Frequency distribution')
    
end; clear indTrial;
