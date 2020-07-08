function eBOSC_example()

% Conduct BOSC + eBOSC analysis

% 180201 | CSD transformation prior to eBOSC calculation

clear all; clc; restoredefaultpath;

%% set paths

pn.root = '/Volumes/LNDG/Projects/StateSwitch/dynamic/data/eeg/rest/';
pn.dataIn   = [pn.root, 'A_preproc/SA_preproc_study/B_data/C_EEG_FT/'];
pn.dataOut  = [pn.root, 'B_analyses/C_eBOSC_CSD/B_data/A_eBOSCout/']; mkdir(pn.dataOut);
pn.tools    = [pn.root, 'B_analyses/C_eBOSC_CSD/T_tools/'];
addpath('/Volumes/LNDG/Projects/StateSwitch/dynamic/data/eeg/rest/B_analyses/A_MSE_CSD/T_tools/fieldtrip-20170904/'); ft_defaults
%addpath([pn.tools, 'fieldtrip-20170904/']); ft_defaults
addpath(genpath([pn.tools, 'eBOSC-master/'])); 

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

%% load subject list

% N = 47 YA, 52 OAs
% 2201 - no rest available; 1213 dropped (weird channel arrangement)

IDs = {'1117';'1118';'1120';'1124';'1126';'1131';'1132';'1135';'1136';'1138';...
    '1144';'1151';'1158';'1160';'1163';'1164';'1167';'1169';'1172';'1173';...
    '1178';'1182';'1215';'1216';'1219';'1221';'1223';'1227';'1228';...
    '1233';'1234';'1237';'1239';'1240';'1243';'1245';'1247';'1250';'1252';...
    '1257';'1261';'1265';'1266';'1268';'1270';'1276';'1281';...
    '2104';'2107';'2108';'2112';'2118';'2120';'2121';'2123';'2125';'2129';...
    '2130';'2131';'2132';'2133';'2134';'2135';'2139';'2140';'2145';'2147';...
    '2149';'2157';'2160';'2202';'2203';'2205';'2206';'2209';'2210';...
    '2211';'2213';'2214';'2215';'2216';'2217';'2219';'2222';'2224';'2226';...
    '2227';'2236';'2237';'2238';'2241';'2244';'2246';'2248';'2250';'2251';...
    '2252';'2258';'2261'};

%% get info about current subject

for indID = 1:numel(IDs)

    results = [];
    
    %% load data

    load([pn.dataIn,  IDs{indID}, '_rest_EEG_Rlm_Fhl_rdSeg_Art_EO.mat'], 'data')

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
    
    rem             = [];
    rem.channel     = {'all','-IOR','-LHEOG','-RHEOG','-A1','-A2'};
    rem.demean      = 'no';
    
    data = ft_preprocessing(rem, data); clear rem;
   
   %% apply CSD transformation
        
    % CSD transform
    csd_cfg = [];
    csd_cfg.elecfile = 'standard_1005.elc';
    csd_cfg.method = 'spline';
    data = ft_scalpcurrentdensity(csd_cfg, data);

    %%%%%%%%%%%
    %% eBOSC %%
    %%%%%%%%%%%
    
    %% oscillation detection 
    
    for e = 1:60 % electrode loop

        % display progress
        display([IDs{indID} ', channel #' num2str(e)])
        tic
        
    %%  TF analysis for whole signal to prepare background fit
    
        B = [];
        for t = 1:length(data.trial) % trial index
            % get data
            tmp_dat = data.trial{t}(e,:);
            % wavelet transform (NOTE: no check to avoid spectral leakage);
            % apply correction factor
            B.trial{t} = BOSC_tf(tmp_dat,cfg.eBOSC.F,cfg.eBOSC.fsample,cfg.eBOSC.wavenumber);
            clear tmp_dat
        end; clear t
           
    %% condition-specific background fit --> power threshold estimation
    
            
            amountTrials = size(B.trial,2);
            
            BG = [];
            for t = 1:length(data.trial)
                % collect TF values; Remove BGpad at beginning and end.
                % Note that BGpad removes some edge segments to avoid 
                % potential edge artifacts, yet the final segment does not
                % exclusively cover the retention period.
                BG = [BG B.trial{t}(:,cfg.eBOSC.BGpad+1:end-cfg.eBOSC.BGpad)];
            end; clear t
            
            %% standard backgrounds
            
            % background power estimation - standard
            [pv_stnd,~] = BOSC_bgfit(cfg.eBOSC.F,BG); 
            mp_stnd = 10.^(polyval(pv_stnd,log10(cfg.eBOSC.F))); 

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

            % BOSC thresholds
            [pt,dt] = BOSC_thresholds(cfg.eBOSC.fsample,cfg.eBOSC.percentile,cfg.eBOSC.ncyc,cfg.eBOSC.F,mp);

            % keep overall background, 1/f fit, and power threshold
            BGinfo.all.(['bg_pow'])(e,:)        = mean(BG(:,cfg.eBOSC.trialPad+1:end-cfg.eBOSC.trialPad),2);
            BGinfo.all.(['bg_log10_pow'])(e,:)  = mean(log10(BG(:,cfg.eBOSC.trialPad+1:end-cfg.eBOSC.trialPad)),2);
            BGinfo.all.(['bg_amp'])(e,:)        = mean(sqrt(BG(:,cfg.eBOSC.trialPad+1:end-cfg.eBOSC.trialPad)),2);
            BGinfo.all.(['pv'])(e,:)            = pv;
            BGinfo.all.(['mp'])(e,:)            = mp;
            BGinfo.all.(['mp_stnd'])(e,:)       = mp_stnd;
            BGinfo.all.(['pt'])(e,:)            = pt;

            % clear variables
            clear pv pv_stnd mp_stnd
            
            for t = 1:size(B.trial,2)

                % initialize variables
                bosc.pepisode{1,t}(e,:)  = zeros(1,size(cfg.eBOSC.F,2));
                bosc.abundance{1,t}(e,:) = zeros(1,size(cfg.eBOSC.F,2));
                bosc.episodes{1,t}{e,1}  = [];

                % wavelet transform
                B_ = B.trial{1,t}; 

            %%  detect rhythms + calculate Pepisode (~standard BOSC)

                time = data.time{1,1}(:,cfg.eBOSC.trialPad+1:end-cfg.eBOSC.trialPad);

                % BOSC oscillation detection
                % WLpadding is removed to avoid edge artifacts during the
                % detection. Note that detectedPad still remains so that there
                % is no problems with too few sample points at the edges to
                % fulfill the numcycles criterion.

                % oscillation detection
                detected = zeros(size(B_(:,cfg.eBOSC.WLpadding+1:end-cfg.eBOSC.WLpadding)));
                for f = 1:length(cfg.eBOSC.F)
                    detected(f,:) = BOSC_detect(B_(f,cfg.eBOSC.WLpadding+1:end-cfg.eBOSC.WLpadding),pt(f),dt(f),cfg.eBOSC.fsample);
                end; clear f

                alphaDetected = zeros(1,size(detected,2));
                alphaDetected(nanmean(detected(cfg.eBOSC.F > 8 & cfg.eBOSC.F < 12,:),1)>0) = 1;
                
                origData = data.trial{1}(e, cfg.eBOSC.WLpadding+1:end-cfg.eBOSC.WLpadding);
                
                results.detectedAlpha(e,:) = alphaDetected;
                results.origData(e,:) = origData;

        end; clear t;
    
    end; clear e; % electrode loop

save([pn.dataOut IDs{indID}, '_MatrixDetected_v1.mat'],'results', 'cfg')

% figure; imagesc(zscore(results.origData,[],2).*results.detectedAlpha) % temporally anticorrelated alpha
% figure; imagesc(results.origData)
% 
% RMat = corrcoef((results.origData.*results.detectedAlpha)');
% figure; imagesc(RMat)
% 
% RMat = corrcoef((results.origData)');
% figure; imagesc(RMat)
% 
% figure; hold on; 
% plot(zscore(results.origData(31,:).*results.detectedAlpha(31,:),[],2));
% plot(zscore(results.origData(33,:).*results.detectedAlpha(33,:),[],2));
% plot(zscore(results.origData(36,:).*results.detectedAlpha(36,:),[],2));

end % function end