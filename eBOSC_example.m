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
cfg.eBOSC.threshold.duration	= repmat(0, 1, numel(cfg.eBOSC.F));         % vector of duration thresholds at each frequency (previously: ncyc)
cfg.eBOSC.threshold.percentile  = .95;                                      % percentile of background fit for power threshold

% episode post-processing
cfg.eBOSC.postproc.use      = 'no';         % Post-processing of rhythmic eBOSC.episodes, i.e., wavelet 'deconvolution' (default = 'no')
cfg.eBOSC.postproc.method   = 'MaxBias';	% Deconvolution method (default = 'MaxBias', FWHM: 'FWHM')
cfg.eBOSC.postproc.edgeOnly = 'yes';        % Deconvolution only at on- and offsets of eBOSC.episodes? (default = 'yes')
cfg.eBOSC.postproc.effSignal= 'PT';         % Amplitude deconvolution on whole signal or signal above power threshold? (default = 'PT')

% general processing settings
cfg.eBOSC.channel = [59,60]; % select channels (default: all)
cfg.eBOSC.trial = [1]; % select trials (default: all)
cfg.eBOSC.trial_background = [1]; % select trials for background (default: all)

%% load data

load([pn.eBOSC,  'util/1160_rest_EEG_Rlm_Fhl_rdSeg_Art_EO.mat'], 'data')

%% concatenate trials for resting state here

data.trial{1} = cat(2,data.trial{:}); data.trial(2:end) = [];
data.time{1} = cat(2,data.time{:}); data.time(2:end) = [];
data = rmfield(data, 'sampleinfo');

%% run eBOSC

[eBOSC, cfg] = eBOSC_wrapper(cfg, data);

%% multiple figures as visual sanity-checks

% Supplementary Figure: plot estimated background + power threshold

indChan = 1; indTrial = 1;
% Note: indChan indexes the requested channel in cfg.eBOSC.channel

figure; hold on;
plot(log10(cfg.eBOSC.F), log10(eBOSC.static.mp(cfg.eBOSC.channel(indChan),:)), 'k--','LineWidth', 1.5); 
plot(log10(cfg.eBOSC.F), log10(eBOSC.static.pt(cfg.eBOSC.channel(indChan),:)), 'k-', 'LineWidth', 1.5)
plot(log10(cfg.eBOSC.F),log10(eBOSC.static.bg_pow(cfg.eBOSC.channel(indChan),:)), 'r-', 'LineWidth', 2)
xlabel('Frequency (log10 Hz)'); ylabel('Power (log 10 a.u.)');
legend({'Aperiodic fit', 'Statistical power threshold', 'Avg. spectrum'}, ...
    'orientation', 'vertical', 'location', 'SouthWest'); legend('boxoff');
set(findall(gcf,'-property','FontSize'),'FontSize',20)
xlim([.3, 1.75])

% Supplementary Plot: plot only rhythmic episodes
h = figure('units','normalized','position',[.1 .1 .6 .3]);
hold on; 
plot(squeeze(eBOSC.origData(indChan, indTrial,:)), 'k');
tmpDetected = squeeze(eBOSC.detectedAlpha(indChan, indTrial,:)); tmpDetected(tmpDetected==0) = NaN;
plot(squeeze(eBOSC.origData(indChan, indTrial,:)).*tmpDetected, 'r');
xlim([7.2, 7.9]*10^4)
xlabel('Time (s)'); ylabel('Amplitude [µV]');
legend({'Original signal'; 'Rhythmic signal'}, ...
    'orientation', 'horizontal', 'location', 'north'); legend('boxoff')
set(findall(gcf,'-property','FontSize'),'FontSize',26)

% Supplementary Plots: different statistics
figure; 
subplot(3,2,1); histogram(eBOSC.episodes.SNRMean); title('SNR distribution')
subplot(3,2,2); histogram(log10(eBOSC.episodes.SNRMean)); title('SNR distribution(log10)')
subplot(3,2,3); histogram(eBOSC.episodes.DurationC); title('Duration distribution')
subplot(3,2,4); histogram(log10(eBOSC.episodes.DurationC)); title('Duration distribution(log10)')
subplot(3,2,5); histogram(eBOSC.episodes.FrequencyMean); title('Frequency distribution')
subplot(3,2,6); hold on; plot(squeeze(eBOSC.pepisode(indChan, indTrial,:))); plot(squeeze(eBOSC.abundance_ep(indChan, indTrial,:))); title('Pepisode, abundance')

% Supplementary Plot: plot rhythmic episodes with indicated onsets (here across the two channels)
idx_alpha = find(eBOSC.episodes.FrequencyMean > 8 & eBOSC.episodes.FrequencyMean <15); % filter for alpha
for indEp = 1:numel(idx_alpha)
    % These are two alternative ways to extract the onset timepoint
    % from the table
    idx_onsetTime(indEp) = find(cfg.tmp.finalTime>= eBOSC.episodes.Onset(idx_alpha(indEp)), 1, 'first');
    idx_onset(indEp) = eBOSC.episodes.ColID{idx_alpha(indEp)}(1);
end
figure; hold on; 
scatter(idx_onset, repmat(100,1,numel(idx_onset)), 'filled')
OnsetLine = zeros(size(squeeze(eBOSC.origData(indChan, indTrial,:))));
OnsetLine(idx_onset) = 100;
plot(OnsetLine, 'g')
plot(squeeze(eBOSC.origData(indChan, indTrial,:)), 'k');
tmpDetected = squeeze(eBOSC.detectedAlpha_ep(indChan, indTrial,:)); tmpDetected(tmpDetected==0) = NaN;
plot(squeeze(eBOSC.origData(indChan, indTrial,:)).*tmpDetected, 'r');
xlim([7.2, 7.9]*10^4)

% General note: the above onset display is not expected to be perfect as an episode with a
% mean frequency within the requested range could have some time points
% falling out of that range.