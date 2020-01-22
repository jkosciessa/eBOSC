%% Plot example of eBOSC episode detection (i.e. post-processing)

% pn.root = '/Volumes/EEG/BOSC_Simulation/17_JQK_RED/A_scripts/K_simulation3c/T_tools/eBOSC-master/'; addpath(genpath(pn.root));
% pn.backgroundData  = '/Volumes/EEG/BOSC_Simulation/17_JQK_RED/A_scripts/K_simulation3c/B_data/'; % INDICATE LOCATION OF SIMULATED BACKGROUND
% pn.out = '/Volumes/EEG/BOSC_Simulation/17_JQK_RED/A_scripts/K_simulation3c/B_data/'; % CHOOSE OUTPUT DIRECTORY
% 
% variant = 'B';
% 
% %% encode run date
% 
% cfg.runDate = date;
% 
% %% simulation parameters
% 
% cfg.simParams.amplitude     = [0 2 4 6 8 12 16 24];                         % simulated signal power
% cfg.simParams.cycles        = [2 4 8 16 32 64 128 200];                     % simulated signal durations [in cycles]; total of 14 seconds
% cfg.simParams.segmentDur    = 20;
% cfg.simParams.time          = [.004:.004:cfg.simParams.segmentDur];         % time vector of complete segment
% cfg.simParams.repetitions   = 50;                                          % amount of repetitions
% 
% %%  eBOSC parameters
% 
% cfg.eBOSC.F                 = 2.^[1:.125:5.25];                             % setup (Whitten et al., 2011), but higher frequency resolution
% cfg.eBOSC.wavenumber        = 3;
% cfg.eBOSC.ncyc              = repmat(3, 1, numel(cfg.eBOSC.F));
% cfg.eBOSC.percentile        = .95;
% cfg.eBOSC.fsample           = 250;
% cfg.eBOSC.WLpadding         = 500;                                          % padding to avoid edge artifacts due to WL [SPs]
% cfg.eBOSC.detectedPad       = 250;                                          % 'shoulder' for BOSC detected matrix to account for duration threshold
% cfg.eBOSC.trialPad          = 750;                                          % complete padding (WL + shoulder)
% cfg.eBOSC.BGpad             = 750;                                          % padding of segments for BG (only avoiding edge artifacts)
% cfg.eBOSC.fres              = 1.2;                                          % cf. Linkenkaer-Hansen, K., et al. (2001). "Long-Range Temporal Correlations and Scaling Behavior in Human Brain Oscillations." The Journal of Neuroscience 21(4): 1370-1377.
% cfg.eBOSC.fstp              = 1;
% cfg.eBOSC.freqRemoval       = 'JQK';
% cfg.eBOSC.BiasCorrection    = 'yes';
% cfg.eBOSC.LowFreqExcludeBG  = 8;
% cfg.eBOSC.HighFreqExcludeBG = 15;
% 
% %%  initialize output matrices
% 
% load([pn.backgroundData 'background.mat'],'bckgrnd_filt')
% 
% amountAmps          = numel(cfg.simParams.amplitude);
% amountCycles        = numel(cfg.simParams.cycles);
% amountTimePoints    = numel(cfg.simParams.time);
% amountFreqs         = numel(cfg.eBOSC.F);
% amountRepetitions   = cfg.simParams.repetitions;
% 
% SignalDetection.Hits    = NaN(amountAmps, amountCycles, amountRepetitions,1);
% SignalDetection.Misses  = NaN(amountAmps, amountCycles, amountRepetitions,1);
% SignalDetection.CAs     = NaN(amountAmps, amountCycles, amountRepetitions,1);
% SignalDetection.FAs     = NaN(amountAmps, amountCycles, amountRepetitions,1);
% 
% abundance_spec_ep   = NaN(amountAmps,amountCycles,amountRepetitions,amountFreqs);
% abundance_ep        = NaN(amountAmps,amountCycles,amountRepetitions,1);
% abundance_spec      = NaN(amountAmps,amountCycles,amountRepetitions, amountFreqs);
% 
% %%  loop conditions
% 
% a = 8;
% c = 3;
%         
% %% create data segments according to power and duration
% for k = 1:amountRepetitions
%     % generate alpha in the middle of the segment
%     rhythmCycles = cfg.simParams.cycles(c);
%     rhythmFreq = 10;
%     rhythmTime(c) = round((rhythmCycles/rhythmFreq),3);
%     % simulate rhythms as symmetrical around the center
%     timeNew = round(rhythmTime(c)/0.004,0);
%     if mod(timeNew,2) ~= 0
%         timeNew = timeNew + 1;
%         rhythmTime(c) = timeNew.*0.004;
%     else rhythmTime(c) = timeNew.*0.004;
%     end; clear timeNew;
%     rhythmTimeVector = [.004:.004:rhythmTime(c)];
%     rhythmIdxVector = (amountTimePoints/2)-(numel(rhythmTimeVector)/2)+1:(amountTimePoints/2)+(numel(rhythmTimeVector)/2);
%     % filter entire signal between 8 and 12 Hz (6th order butterworth) (not locally on alpha)
%     [tmp_b,tmp_a] = butter(6, [8, 12]/(250/2), 'bandpass');
%     tmp_bpsignal = filter(tmp_b,tmp_a,bckgrnd_filt(k,:));
%     %VarBG = var(tmp_bpsignal(rhythmIdxVector));
%     VarBG = var(tmp_bpsignal);
%     % scale rhythm by noise background power
%     targetPower = cfg.simParams.amplitude(a)*VarBG;
%     amplitudeFromRMS = (sqrt(targetPower)*sqrt(2));
%     simulatedRhythm = sin(rhythmTimeVector*2*pi*rhythmFreq)*amplitudeFromRMS;
%     simulatedRhythm_complete = zeros(1,amountTimePoints);
%     simulatedRhythm_complete(1,rhythmIdxVector) = simulatedRhythm;
%     % effective segment = BG + rhythm
%     signal = bckgrnd_filt(k,:) + simulatedRhythm_complete;
%     %  wavelet transform to inspect 'simulated SNR'
%     B(k,:,:) = BOSC_tf(signal,cfg.eBOSC.F,cfg.eBOSC.fsample,cfg.eBOSC.wavenumber);
%     BGonly_FT(k,:,:) = BOSC_tf(bckgrnd_filt(k,:),cfg.eBOSC.F,cfg.eBOSC.fsample,cfg.eBOSC.wavenumber);
%     freq = 19;
%     FT_SNR_local(a,c,k) = mean(B(k,freq,rhythmIdxVector))./mean(BGonly_FT(k,freq,rhythmIdxVector));
%     FT_SNR_total(a,c,k) = mean(B(k,19,rhythmIdxVector))./mean(BGonly_FT(k,19,:));
% end % k repetitions
% 
% %%  run extended BOSC BG estimation
% 
% % background is estimated across all simulated trials
% BG = [];
% for k = 1:amountRepetitions
%     BG = cat(3,BG, B(k,:,:));
% end
% BG = squeeze(BG);
% 
% % background power estimation - robust
% % find peak between 8-15 Hz, get wavelet extension in frequency domain, remove
% % points within this range from the estimation; compute robust
% % regression
% 
% freqInd1 = find(cfg.eBOSC.F >= cfg.eBOSC.LowFreqExcludeBG, 1, 'first');
% freqInd2 = find(cfg.eBOSC.F <= cfg.eBOSC.HighFreqExcludeBG, 1, 'last');
% 
% [~, indPos] = max(mean(BG(freqInd1:freqInd2,:),2));
% indPos = freqInd1+indPos;
% 
% LowFreq = cfg.eBOSC.F(indPos)-(((2/cfg.eBOSC.wavenumber)*cfg.eBOSC.F(indPos))/2);
% UpFreq = cfg.eBOSC.F(indPos)+(((2/cfg.eBOSC.wavenumber)*cfg.eBOSC.F(indPos))/2);
% 
% freqIndLow = find(cfg.eBOSC.F >= LowFreq, 1, 'first');
% freqIndHigh = find(cfg.eBOSC.F <= UpFreq, 1, 'last');
% 
% [pv,~] = eBOSC_bgfit_robust(cfg.eBOSC.F([1:freqIndLow-1 freqIndHigh+1:end]),BG([1:freqIndLow-1 freqIndHigh+1:end], :));
% mp = 10.^(polyval(pv,log10(cfg.eBOSC.F)));
% 
% % thresholds
% [pt,dt] = BOSC_thresholds(cfg.eBOSC.fsample,cfg.eBOSC.percentile,cfg.eBOSC.ncyc,cfg.eBOSC.F,mp);
% 
% % keep overall background, 1/f fit, and power threshold
% BGinfo.all.bg_pow(a,c,k,:) = mean(BG(:,cfg.eBOSC.trialPad+1:end-cfg.eBOSC.trialPad),2);
% BGinfo.all.bg_log10_pow(a,c,k,:) = mean(log10(BG(:,cfg.eBOSC.trialPad+1:end-cfg.eBOSC.trialPad)),2);
% BGinfo.all.bg_amp(a,c,k,:) = mean(sqrt(BG(:,cfg.eBOSC.trialPad+1:end-cfg.eBOSC.trialPad)),2);
% BGinfo.all.pv(a,c,k,:) = pv;
% BGinfo.all.mp(a,c,k,:) = mp;
% BGinfo.all.pt(a,c,k,:)= pt;
% 
% %% Plot example of eBOSC processing (gather Figure data)
% 
% k = 1;
% 
% % oscillation detection
% detected = zeros(35,amountTimePoints-2*cfg.eBOSC.WLpadding);
% for f = 1:length(cfg.eBOSC.F)
%     detected(f,:) = BOSC_detect(squeeze(B(k,f,cfg.eBOSC.WLpadding+1:end-cfg.eBOSC.WLpadding))',pt(f),dt(f),cfg.eBOSC.fsample);
% end; clear f
% 
% signal = bckgrnd_filt(k,:) + simulatedRhythm_complete;
% 
% Figure2.A.signal = signal;
% Figure2.A.time = [1850+cfg.eBOSC.WLpadding 2200+cfg.eBOSC.WLpadding];
% Figure2.B.signal = simulatedRhythm_complete;
% Figure2.B.time = [1850+cfg.eBOSC.WLpadding 2200+cfg.eBOSC.WLpadding];
% Figure2.C.detected = detected;
% Figure2.C.freq = cfg.eBOSC.F;
% 
% cfg.eBOSC.npnts = size(detected,2);
% cfg.eBOSC.pt  = pt;
% cfg.eBOSC.BiasVersion = '171023_MaxBias_PT';
% cfg.eBOSC.BiasCorrection = 'yes';
% cfg.eBOSC.method = 'MaxBias';
% cfg.eBOSC.edgeOnly = 'yes';
% cfg.eBOSC.effSignal = 'PT';
% 
% addpath('/Volumes/EEG/BOSC_Simulation/17_JQK_RED/A_scripts/L_ExampleEpisodeCreation/A_scripts/');
% 
% episodes = [];
% [detected1,episodes,Figure2] = eBOSC_createEpisodes_plotVersion(squeeze(B(k,:,cfg.eBOSC.WLpadding+1:end-cfg.eBOSC.WLpadding)),detected, cfg, h, Figure2);
% 
% %% save Figure data
% 
% save('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/F2.mat', 'Figure2')

%% load Figure data

load('/Users/kosciessa/Desktop/eBOSC/figureData/S1.mat', 'Figure2')

%% plot Figure

h = figure('units','normalized','position',[.1 .1 .7 .7]);
subplot(3,2,1); plot(Figure2.A.signal, 'k', 'LineWidth', 2); xlim(Figure2.A.time); title('Time domain signal')
    set(gca, 'XTick', [2400 2600]); set(gca, 'XTickLabel', {'Alpha ON'; 'Alpha OFF'}); ylim([-7 7]); ylabel({'Time-domain','amplitude (a.u.)'});
subplot(3,2,2); plot(Figure2.B.signal, 'k', 'LineWidth', 2); xlim(Figure2.B.time); title('Alpha rhythm only')
    set(gca, 'XTick', [2400 2600]); set(gca, 'XTickLabel', {'Alpha ON'; 'Alpha OFF'}); ylabel({'Time-domain','amplitude (a.u.)'});
subplot(3,2,3); imagesc(Figure2.C.detected); xlim([1850 2200]); title('BOSC detected matrix'); xlabel('Time'); ylabel('Frequency [Hz]')
    set(gca, 'YTick', [10 20 30]); set(gca, 'YTickLabel', round(Figure2.C.freq([10 20 30])));
    set(gca, 'XTick', [1900 2100]); set(gca, 'XTickLabel', {'Alpha ON'; 'Alpha OFF'});
subplot(3,2,4); imagesc(Figure2.D.detected);xlim([1850 2200])
    title('eBOSC sparse continuous matrix')
    ylabel('Frequency [Hz]')
    set(gca, 'YTick', [10 20 30]);
    set(gca, 'YTickLabel', round(Figure2.C.freq([10 20 30])));
    set(gca, 'XTick', [1900 2100]); set(gca, 'XTickLabel', {'Alpha ON'; 'Alpha OFF'});
subplot(3,2,5); hold on;
    plot(Figure2.E.time, Figure2.E.traceBlack, 'k'); hold on; plot(Figure2.E.time, Figure2.E.traceRed, 'r');
    xlim([1850 2200])
    title('Temporal Convolution Correction')
    ylabel({'Frequency-domain','rhythmic amplitude [a.u.]'})
    set(gca, 'XTick', [1900 2100]); set(gca, 'XTickLabel', {'Alpha ON'; 'Alpha OFF'});
subplot(3,2,6); hold on;
    plot(Figure2.F.traceRedX, Figure2.F.traceRedY, 'r', 'LineWidth', 10)
    plot(Figure2.F.traceBlackX, Figure2.F.traceBlackY, 'k', 'LineWidth', 10)
    xlim([1850 2200])
    title('Final Episode Output')
    ylabel({'Frequency-domain',' rhythmic amplitude [a.u.]'})
    set(gca, 'XTick', [1900 2100]); set(gca, 'XTickLabel', {'Alpha ON'; 'Alpha OFF'});
set(findall(gcf,'-property','FontSize'),'FontSize',18)

pn.plotFolder = '/Volumes/fb-lip/BOSC_SternRest/X_documentation/B_2018_Manuscript/F_Figures/';
figureName = 'F2';

saveas(h, [pn.plotFolder, figureName], 'fig');
saveas(h, [pn.plotFolder, figureName], 'epsc');
saveas(h, [pn.plotFolder, figureName], 'png');
saveas(h, [pn.plotFolder, figureName], 'pdf');