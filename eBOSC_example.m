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

%%  BOSC parameters

cfg.eBOSC.F                 = 2.^[1:.125:6];                            % frequency sampling (~Whitten et al., 2011), but higher frequency resolution
cfg.eBOSC.wavenumber        = 6;                                        % wavelet family parameter (time-frequency tradeoff) [recommended: ~6]
cfg.eBOSC.ncyc              = repmat(3, 1, numel(cfg.eBOSC.F));         % vector of duration thresholds at each frequency
cfg.eBOSC.percentile        = .95;                                      % percentile of background fit for power threshold
cfg.eBOSC.fsample           = 500;                                      % current sampling frequency of EEG data
cfg.eBOSC.WLpadding         = 500;                                      % padding to avoid edge artifacts due to WL [SPs]
cfg.eBOSC.detectedPad       = 250;                                      % 'shoulder' for BOSC detected matrix to account for duration threshold
cfg.eBOSC.trialPad          = 750;                                      % complete padding (WL + shoulder)
cfg.eBOSC.BGpad             = 750;                                      % padding of segments for BG (only avoiding edge artifacts)
cfg.eBOSC.fres              = 1.2;                                      % cf. Linkenkaer-Hansen, K., et al. (2001). "Long-Range Temporal Correlations and Scaling Behavior in Human Brain Oscillations." The Journal of Neuroscience 21(4): 1370-1377.
cfg.eBOSC.fstp              = 1;
cfg.eBOSC.freqRemoval       = 'JQK';
cfg.eBOSC.BiasCorrection    = 'yes';                                    % use temporal correction for impact of wavelet?
cfg.eBOSC.method            = 'MaxBias';
cfg.eBOSC.edgeOnly          = 'no';
cfg.eBOSC.effSignal         = 'PT';
cfg.eBOSC.LowFreqExcludeBG  = 8;                                        % lower bound of bandpass to be excluded prior to background fit
cfg.eBOSC.HighFreqExcludeBG = 15;                                       % higher bound of bandpass to be excluded prior to background fit
cfg.eBOSC.waveseg           = [0 9];                                    % include +-3s around stim processing; 3 seconds will be cut at each end during detection --> 3 to 6 (stim only)

% -----------------------------
% TO DO: rename BiasCorrection to postproc, turn into structure
% -----------------------------

% cfg.eBOSC.postproc.use               = 'yes'; (default = 'no')
% cfg.eBOSC.postproc.method            = 'MaxBias'; (default = 'MaxBias', FWHM: 'FH')
% cfg.eBOSC.postproc.edgeOnly          = 'no'; (default = 'yes')
% cfg.eBOSC.postproc.effSignal         = 'PT'; (default = 'PT')

eBOSC = [];

%% load data

load([pn.eBOSC,  'util/1160_rest_EEG_Rlm_Fhl_rdSeg_Art_EC.mat'], 'data')

% concatenate trials: check whether this is reasonable, as it introduces
% hard cuts

data.trial{1} = cat(2,data.trial{:}); data.trial(2:end) = [];
data.time{1} = cat(2,data.time{:}); data.time(2:end) = [];
data = rmfield(data, 'sampleinfo');

eBOSC.inputTime = data.time{1,1};
eBOSC.detectedTime = eBOSC.inputTime(cfg.eBOSC.WLpadding+1:end-cfg.eBOSC.WLpadding);
eBOSC.finalTime = eBOSC.inputTime(cfg.eBOSC.trialPad+1:end-cfg.eBOSC.trialPad);

%% ---- select a channel here

e = 60;

% -----------------------------
% TO DO: include fieldtrip-style automatic detection of dimord
% -----------------------------

%% oscillation detection 

% display progress
display(['channel #' num2str(e)])
tic

%%  TF analysis for whole signal to prepare background fit

info.Ntrial = length(data.trial);

TFR = [];
for indTrial = 1:info.Ntrial
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
for indTrial = 1:info.Ntrial
    % remove BGpad at beginning and end to avoid edge artifacts
    BG = [BG TFR.trial{indTrial}(:,cfg.eBOSC.BGpad+1:end-cfg.eBOSC.BGpad)];
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
BGinfo.all.(['bg_pow'])(e,:)        = mean(BG(:,cfg.eBOSC.trialPad+1:end-cfg.eBOSC.trialPad),2);
% log10-transformed wavelet power spectrum (NOT only background)
BGinfo.all.(['bg_log10_pow'])(e,:)  = mean(log10(BG(:,cfg.eBOSC.trialPad+1:end-cfg.eBOSC.trialPad)),2);
% intercept and slope parameters of the robust linear 1/f fit (log-log)
BGinfo.all.(['pv'])(e,:)            = pv;
% linear background power at each estimated frequency
BGinfo.all.(['mp'])(e,:)            = mp;
% statistical power threshold
BGinfo.all.(['pt'])(e,:)            = pt;

%% TO DO: implement IRASA for background estimation

%% use statically-defined threshold for single-trial detection

for indTrial = 1:info.Ntrial

    % initialize variables
    eBOSC.pepisode{1,indTrial}(e,:)  = zeros(1,size(cfg.eBOSC.F,2));
    eBOSC.abundance{1,indTrial}(e,:) = zeros(1,size(cfg.eBOSC.F,2));
    eBOSC.episodes{1,indTrial}{e,1}  = [];

    % get wavelet transform for single trial
    % WLpadding is removed to avoid edge artifacts during the
    % detection. Note that detectedPad still remains so that there
    % is no problems with too few sample points at the edges to
    % fulfill the numcycles criterion.
    TFR_ = TFR.trial{1,indTrial}(:,cfg.eBOSC.WLpadding+1:end-cfg.eBOSC.WLpadding);

%%  Detect rhythms + calculate Pepisode (~standard BOSC)

    % The next section applies both the power and the optional duration
    % threshold to detect individual rhythmic segments in the continuous signals.

    detected = zeros(size(TFR_));
    for f = 1:length(cfg.eBOSC.F)
        detected(f,:) = BOSC_detect(TFR_(f,:),pt(f),dt(f),cfg.eBOSC.fsample);
    end; clear f

    % encode detected alpha signals (optional)
    alphaDetected = zeros(1,size(detected,2));
    alphaDetected(nanmean(detected(cfg.eBOSC.F > 8 & cfg.eBOSC.F < 12,:),1)>0) = 1;
    eBOSC.detectedAlpha(e,:) = alphaDetected;

    % encode original signals (optional)
    origData = data.trial{1}(e, cfg.eBOSC.WLpadding+1:end-cfg.eBOSC.WLpadding);
    eBOSC.origData(e,:) = origData;

    % Supplementary Plot: plot only rhythmic episodes
    % for a great example see e.g., Trial 1

    figure; hold on; 
    plot(eBOSC.origData(e,:));
    tmpDetected = eBOSC.detectedAlpha; tmpDetected(tmpDetected==0) = NaN;
    plot(eBOSC.origData(e,:).*tmpDetected(e,:));
    xlim([7.2, 7.9]*10^4)

    %% create list of rhythmic episodes

    % add size of the complete detected matrix and power threshold
    cfg.eBOSC.npnts = size(detected,2);
    cfg.eBOSC.pt  = pt;

    % -----------------------------
    % TO DO: improve episode output 
    % -----------------------------

    episodes = [];
    [detected1,episodes] = eBOSC_createEpisodes(squeeze(TFR_),detected, cfg);

    % encode abundance of episodes (optional)
    eBOSC.abundance_ep(a,c,k,:) = mean(detected1(:,cfg.eBOSC.detectedPad+1:end-cfg.eBOSC.detectedPad),2);

    % Supplementary Plot: original detected vs. sparse episode power
    figure; 
    subplot(121); imagesc(detected.*squeeze(TFR_));
    subplot(122); imagesc(detected1.*squeeze(TFR_));


end; clear indTrial;
