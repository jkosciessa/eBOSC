% Example script for eBOSC analysis

%% initial clean-up

clear all; clc; restoredefaultpath;

%% add eBOSC toolbox with subdirectories

% automatically get script location from editor handle
tmp = matlab.desktop.editor.getActive;
pn.eBOSC = [fileparts(tmp.Filename), '/']; clear tmp;

cd(pn.eBOSC)
addpath(pn.eBOSC);
addpath([pn.eBOSC, 'external/BOSC/']);
addpath([pn.eBOSC, 'external/fieldtrip/']); ft_defaults;
addpath([pn.eBOSC, 'external/NoiseTools/']);

%%  eBOSC parameters

% general setup
cfg.eBOSC.F                 = 2.^[1:.125:6];                            % frequency sampling (~Whitten et al., 2011), but higher frequency resolution
cfg.eBOSC.wavenumber        = 6;                                        % wavelet family parameter (time-frequency tradeoff) [recommended: ~6]
cfg.eBOSC.fsample           = 500;                                      % current sampling frequency of EEG data

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
cfg.eBOSC.ncyc              = repmat(3, 1, numel(cfg.eBOSC.F));         % vector of duration thresholds at each frequency
cfg.eBOSC.percentile        = .95;                                      % percentile of background fit for power threshold
cfg.eBOSC.LowFreqExcludeBG  = 8;                                        % lower bound of bandpass to be excluded prior to background fit
cfg.eBOSC.HighFreqExcludeBG = 15;                                       % higher bound of bandpass to be excluded prior to background fit

% episode creation
cfg.eBOSC.fstp              = 1;

% episode post-processing
cfg.eBOSC.postproc.use      = 'yes';                                    % Post-processing of rhythmic eBOSC.episodes, i.e., wavelet 'deconvolution' (default = 'no')
cfg.eBOSC.postproc.method   = 'MaxBias';                                % Deconvolution method (default = 'MaxBias', FWHM: 'FWHM')
cfg.eBOSC.postproc.edgeOnly = 'no';                                     % Deconvolution only at on- and offsets of eBOSC.episodes? (default = 'yes')
cfg.eBOSC.postproc.effSignal= 'PT';                                     % Amplitude deconvolution on whole signal or signal above power threshold? (default = 'PT')

%% load data

load([pn.eBOSC,  'util/1160_rest_EEG_Rlm_Fhl_rdSeg_Art_EC.mat'], 'data')

%% concatenate trials for resting state here
% TO DO: check whether this is reasonable, as it introduces hard cuts

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

% average power estimates across periods of interest
BG = [];
for indTrial = 1:eBOSC.Ntrial
    % remove BGpad at beginning and end to avoid edge artifacts
    BG = [BG TFR.trial{indTrial}(:,cfg.eBOSC.pad.background_sample+1:end-cfg.eBOSC.pad.background_sample)];
end; clear indTrial

% find peak between 8-15 Hz
freqInd1 = find(cfg.eBOSC.F >= cfg.eBOSC.LowFreqExcludeBG, 1, 'first');
freqInd2 = find(cfg.eBOSC.F <= cfg.eBOSC.HighFreqExcludeBG, 1, 'last');
[~, indPos] = max(mean(BG(freqInd1:freqInd2,:),2));
indPos = freqInd1+indPos;
% approximate wavelet extension in frequency domain
LowFreq = cfg.eBOSC.F(indPos)-(((2/cfg.eBOSC.wavenumber)*cfg.eBOSC.F(indPos))/2);
UpFreq = cfg.eBOSC.F(indPos)+(((2/cfg.eBOSC.wavenumber)*cfg.eBOSC.F(indPos))/2);
% remove power estimates within the following range from the fit
freqIndLow = find(cfg.eBOSC.F >= LowFreq, 1, 'first');
freqIndHigh = find(cfg.eBOSC.F <= UpFreq, 1, 'last');
% perform the robust linear fit, only including putatively aperiodic components (i.e., peak exclusion)
[pv,~] = eBOSC_bgfit_robust(cfg.eBOSC.F([1:freqIndLow-1 freqIndHigh+1:end]),BG([1:freqIndLow-1 freqIndHigh+1:end], :));
mp = 10.^(polyval(pv,log10(cfg.eBOSC.F))); 

% compute eBOSC power (pt) and duration (dt) thresholds: 
[pt,dt] = BOSC_thresholds(cfg.eBOSC.fsample,cfg.eBOSC.percentile,cfg.eBOSC.ncyc,cfg.eBOSC.F,mp);

% save multiple time-invariant estimates that could be of interest:
% overall wavelet power spectrum (NOT only background)
eBOSC.static.bg_pow(e,:)        = mean(BG(:,cfg.eBOSC.pad.total_sample+1:end-cfg.eBOSC.pad.total_sample),2);
% log10-transformed wavelet power spectrum (NOT only background)
eBOSC.static.bg_log10_pow(e,:)  = mean(log10(BG(:,cfg.eBOSC.pad.total_sample+1:end-cfg.eBOSC.pad.total_sample)),2);
% intercept and slope parameters of the robust linear 1/f fit (log-log)
eBOSC.static.pv(e,:)            = pv;
% linear background power at each estimated frequency
eBOSC.static.mp(e,:)            = mp;
% statistical power threshold
eBOSC.static.pt(e,:)            = pt;

%% TO DO: implement IRASA for background estimation

%% use statically-defined threshold for single-trial detection

for indTrial = 1:eBOSC.Ntrial

    % initialize variables
    eBOSC.pepisode{1,indTrial}(e,:)  = zeros(1,size(cfg.eBOSC.F,2));
    eBOSC.abundance{1,indTrial}(e,:) = zeros(1,size(cfg.eBOSC.F,2));
    eBOSC.eBOSC.episodes{1,indTrial}{e,1}  = [];
        
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

    % encode eBOSC.detected alpha signals (optional)
    alphaDetected = zeros(1,size(eBOSC.detected,2));
    alphaDetected(nanmean(eBOSC.detected(cfg.eBOSC.F > 8 & cfg.eBOSC.F < 12,:),1)>0) = 1;
    eBOSC.detectedAlpha(e,:) = alphaDetected;

    % encode original signals (optional)
    origData = data.trial{1}(e, cfg.eBOSC.pad.tfr_sample+1:end-cfg.eBOSC.pad.tfr_sample);
    eBOSC.origData(e,:) = origData;

    % Supplementary Plot: plot only rhythmic eBOSC.episodes
    figure; hold on; 
    plot(eBOSC.origData(e,:));
    tmpDetected = eBOSC.detectedAlpha; tmpDetected(tmpDetected==0) = NaN;
    plot(eBOSC.origData(e,:).*tmpDetected(e,:));
    xlim([7.2, 7.9]*10^4)

    %% create list of rhythmic eBOSC.episodes

    eBOSC.episodes = [];
    [eBOSC.detected1,eBOSC.episodes] = eBOSC_episode_create(squeeze(TFR_),eBOSC, cfg);

    % encode abundance of eBOSC.episodes (optional)
    eBOSC.abundance_ep(a,c,k,:) = mean(eBOSC.detected1(:,cfg.eBOSC.pad.detection_sample+1:end-cfg.eBOSC.pad.detection_sample),2);

    % Supplementary Plot: original eBOSC.detected vs. sparse episode power
    figure; 
    subplot(121); imagesc(eBOSC.detected.*squeeze(TFR_));
    subplot(122); imagesc(eBOSC.detected1.*squeeze(TFR_));
    
    % remove padding for detection (already done for eBOSC.episodes)
    eBOSC.detected1 = eBOSC.detected1(:,cfg.eBOSC.pad.detection_sample+1:end-cfg.eBOSC.pad.detection_sample);

    % Supplementary Plots:

    figure; histogram(log10(eBOSC.episodes.SNRMean))
    figure; histogram(log10(eBOSC.episodes.DurationC))
    
end; clear indTrial;
