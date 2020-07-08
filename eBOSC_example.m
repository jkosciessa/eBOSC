function eBOSC_example()

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
cfg.eBOSC.BiasVersion       = '170630_MaxBias_PT';
cfg.eBOSC.BiasCorrection    = 'yes';                                    % use temporal correction for impact of wavelet?
cfg.eBOSC.method            = 'MaxBias';
cfg.eBOSC.edgeOnly          = 'no';
cfg.eBOSC.effSignal         = 'PT';
cfg.eBOSC.LowFreqExcludeBG  = 8;                                        % lower bound of bandpass to be excluded prior to background fit
cfg.eBOSC.HighFreqExcludeBG = 15;                                       % higher bound of bandpass to be excluded prior to background fit
cfg.eBOSC.waveseg           = [0 9];                                    % include +-3s around stim processing; 3 seconds will be cut at each end during detection --> 3 to 6 (stim only)

results = [];

%% load data

load([pn.eBOSC,  'util/1160_rest_EEG_Rlm_Fhl_rdSeg_Art_EC.mat'], 'data')

% concatenate trials: check whether this is reasonable, as it introduces
% hard cuts

data.trial{1} = cat(2,data.trial{:});
data.trial(2:end) = [];
data.time{1} = cat(2,data.time{:});
data.time(2:end) = [];
data = rmfield(data, 'sampleinfo');

bosc.inputTime = data.time{1,1};
bosc.detectedTime = bosc.inputTime(cfg.eBOSC.WLpadding+1:end-cfg.eBOSC.WLpadding);
bosc.finalTime = bosc.inputTime(cfg.eBOSC.trialPad+1:end-cfg.eBOSC.trialPad);

%% remove extraneous channels
%     
%     rem             = [];
%     rem.channel     = {'all','-IOR','-LHEOG','-RHEOG','-A1','-A2'};
%     rem.demean      = 'no';
%     
%     data = ft_preprocessing(rem, data); clear rem;
%    
%    %% apply CSD transformation
%         
%     % CSD transform
%     csd_cfg = [];
%     csd_cfg.elecfile = 'standard_1005.elc';
%     csd_cfg.method = 'spline';
%     data = ft_scalpcurrentdensity(csd_cfg, data);

%% oscillation detection 

for e = 60 % electrode loop

    % display progress
    display(['channel #' num2str(e)])
    tic

%%  TF analysis for whole signal to prepare background fit

    TFR = [];
    for indTrial = 1:length(data.trial) % trial index
        % get data
        tmp_dat = data.trial{indTrial}(e,:);
        % wavelet transform (NOTE: no check to avoid spectral leakage);
        % apply correction factor
        TFR.trial{indTrial} = BOSC_tf(tmp_dat,cfg.eBOSC.F,cfg.eBOSC.fsample,cfg.eBOSC.wavenumber);
        clear tmp_dat
    end; clear indTrial

%% condition-specific background fit --> power threshold estimation

    info.Ntrial = size(TFR.trial,2);

    BG = [];
    for indTrial = 1:info.Ntrial
        % collect TF values; Remove BGpad at beginning and end.
        % Note that BGpad removes some edge segments to avoid 
        % potential edge artifacts, yet the final segment does not
        % exclusively cover the retention period.
        BG = [BG TFR.trial{indTrial}(:,cfg.eBOSC.BGpad+1:end-cfg.eBOSC.BGpad)];
    end; clear indTrial

    %% eBOSC background: robust fit (2020 NeuroImage paper)

    % robust background power estimation
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

    % compute BOSC thresholds: power threshold (pt) and duration
    % threshold (dt)
    [pt,dt] = BOSC_thresholds(cfg.eBOSC.fsample,cfg.eBOSC.percentile,cfg.eBOSC.ncyc,cfg.eBOSC.F,mp);

    % Here, we save multiple variables that could be of interest later:
    % the overall wavelet power spectrum (NOT only background),

    % overall wavelet power spectrum (NOT only background)
    BGinfo.all.(['bg_pow'])(e,:)        = mean(BG(:,cfg.eBOSC.trialPad+1:end-cfg.eBOSC.trialPad),2);

    % log10-transformed wavelet power spectrum (NOT only background)
    BGinfo.all.(['bg_log10_pow'])(e,:)  = mean(log10(BG(:,cfg.eBOSC.trialPad+1:end-cfg.eBOSC.trialPad)),2);
    BGinfo.all.(['bg_amp'])(e,:)        = mean(sqrt(BG(:,cfg.eBOSC.trialPad+1:end-cfg.eBOSC.trialPad)),2);

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
        bosc.pepisode{1,indTrial}(e,:)  = zeros(1,size(cfg.eBOSC.F,2));
        bosc.abundance{1,indTrial}(e,:) = zeros(1,size(cfg.eBOSC.F,2));
        bosc.episodes{1,indTrial}{e,1}  = [];

        % wavelet transform
        TFR_ = TFR.trial{1,indTrial};

    %%  Detect rhythms + calculate Pepisode (~standard BOSC)

        % WLpadding is removed to avoid edge artifacts during the
        % detection. Note that detectedPad still remains so that there
        % is no problems with too few sample points at the edges to
        % fulfill the numcycles criterion.

        % The next section applies both the power and the optional duration
        % threshold to detect individual rhythmic segments in the continuous signals.
        
        detected = zeros(size(TFR_(:,cfg.eBOSC.WLpadding+1:end-cfg.eBOSC.WLpadding)));
        for f = 1:length(cfg.eBOSC.F)
            detected(f,:) = BOSC_detect(TFR_(f,cfg.eBOSC.WLpadding+1:end-cfg.eBOSC.WLpadding),pt(f),dt(f),cfg.eBOSC.fsample);
        end; clear f

        alphaDetected = zeros(1,size(detected,2));
        alphaDetected(nanmean(detected(cfg.eBOSC.F > 8 & cfg.eBOSC.F < 12,:),1)>0) = 1;

        origData = data.trial{1}(e, cfg.eBOSC.WLpadding+1:end-cfg.eBOSC.WLpadding);

        results.detectedAlpha(e,:) = alphaDetected;
        results.origData(e,:) = origData;

        % Supplementary Plot: plot only rhythmic episodes
        % for a great example see e.g., Trial 1
        
        figure; hold on; 
        plot(results.origData(e,:));
        tmpDetected = results.detectedAlpha; tmpDetected(tmpDetected==0) = NaN;
        plot(results.origData(e,:).*tmpDetected(e,:));
        xlim([7.2, 7.9]*10^4)
        
        %% create list of rhythmic episodes
        
        % check whether the next two entries are indeed necessary
        cfg.eBOSC.npnts = size(detected,2);
        cfg.eBOSC.pt  = pt;
        
        episodes = [];
        [detected1,episodes] = eBOSC_createEpisodes(squeeze(TFR_(:,cfg.eBOSC.WLpadding+1:end-cfg.eBOSC.WLpadding)),detected, cfg);

        abundance_spec_ep(a,c,k,:) = mean(detected1(:,cfg.eBOSC.detectedPad+1:end-cfg.eBOSC.detectedPad),2);

        segLength = size(detected1(:,cfg.eBOSC.detectedPad+1:end-cfg.eBOSC.detectedPad),2);

        % extra: plot detected vs detected1 power
        %figure; subplot(121); imagesc(detected.*squeeze(B(k,:,cfg.eBOSC.WLpadding+1:end-cfg.eBOSC.WLpadding))); subplot(122); imagesc(detected1.*squeeze(B(k,:,cfg.eBOSC.WLpadding+1:end-cfg.eBOSC.WLpadding))); 

        % remove episodes and part of episodes that fall into 'shoulder'
        [episodes] = eBOSC_episode_rm_shoulder(cfg,detected1,episodes);
        
    end; clear indTrial;

end; clear e; % electrode loop

