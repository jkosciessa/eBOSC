function [eBOSC, cfg] = eBOSC_wrapper(cfg, data)
% Main eBOSC wrapper function. Executes eBOSC subfunctions.
%
% Inputs: 
%           cfg | config structure with cfg.eBOSC field:
%                     cfg.eBOSC.F                     | frequency sampling
%                     cfg.eBOSC.wavenumber            | wavelet family parameter (time-frequency tradeoff)
%                     cfg.eBOSC.fsample               | current sampling frequency of EEG data
%                     cfg.eBOSC.pad.tfr_s             | padding following wavelet transform to avoid edge artifacts in seconds (bi-lateral)
%                     cfg.eBOSC.pad.detection_s       | padding following rhythm detection in seconds (bi-lateral); 'shoulder' for BOSC eBOSC.detected matrix to account for duration threshold
%                     cfg.eBOSC.pad.total_s           | complete padding (WL + shoulder)
%                     cfg.eBOSC.pad.background_s      | padding of segments for BG (only avoiding edge artifacts)
%                     cfg.eBOSC.threshold.excludePeak | lower and upper bound of frequencies to be excluded during background fit (Hz) (previously: LowFreqExcludeBG HighFreqExcludeBG)
%                     cfg.eBOSC.threshold.duration    | vector of duration thresholds at each frequency (previously: ncyc)
%                     cfg.eBOSC.threshold.percentile  | percentile of background fit for power threshold
%                     cfg.eBOSC.postproc.use          | Post-processing of rhythmic eBOSC.episodes, i.e., wavelet 'deconvolution' (default = 'no')
%                     cfg.eBOSC.postproc.method       | Deconvolution method (default = 'MaxBias', FWHM: 'FWHM')
%                     cfg.eBOSC.postproc.edgeOnly     | Deconvolution only at on- and offsets of eBOSC.episodes? (default = 'yes')
%                     cfg.eBOSC.postproc.effSignal	  | Power deconvolution on whole signal or signal above power threshold? (default = 'PT')
%                     cfg.eBOSC.channel               | Subset of channels? (default: [] = all)
%                     cfg.eBOSC.trial                 | Subset of trials? (default: [] = all)
%                     cfg.eBOSC.trial_background      | Subset of trials for background? (default: [] = all)
%           data | input time series data in FieldTrip format with:
%                | .trial field: {trial}(channel x time)
%                | .time field: {trial}(channel x time)
%                | .label field: {channelName}
%
% Outputs: 
%           eBOSC | main eBOSC output structure
%               eBOSC.episodes | table of individual rhythmic episodes (see eBOSC_episode_create)
%               eBOSC.detected | binary matrix of detected time-frequency points (prior to episode creation)
%               eBOSC.pepisode | temporal average of detected rhythms (prior to episode creation)
%               eBOSC.detected_ep | binary matrix of detected time-frequency points (following episode creation)
%               eBOSC.abundance_ep | temporal average of detected rhythms (following episode creation)
%           cfg | config structure

    eBOSC = [];
    
    % set some defaults for included channels and trials, if not specified
    if isempty(cfg.eBOSC.channel)
        cfg.eBOSC.channel = 1:numel(data.label);
    end
    if isempty(cfg.eBOSC.trial)
        cfg.eBOSC.trial = 1:numel(data.trial);
    end
    if isempty(cfg.eBOSC.trial_background)
        cfg.eBOSC.trial_background = 1:numel(data.trial);
    end
    
    % calculate the sample points for paddding
    cfg.eBOSC.pad.tfr_sample = cfg.eBOSC.pad.tfr_s.*cfg.eBOSC.fsample;                          % automatic sample point calculation
    cfg.eBOSC.pad.detection_sample = cfg.eBOSC.pad.detection_s.*cfg.eBOSC.fsample;              % automatic sample point calculation
    cfg.eBOSC.pad.total_s = cfg.eBOSC.pad.tfr_s + cfg.eBOSC.pad.detection_s;                    % complete padding (WL + shoulder)
    cfg.eBOSC.pad.total_sample = cfg.eBOSC.pad.tfr_sample + cfg.eBOSC.pad.detection_sample;
    cfg.eBOSC.pad.background_sample = cfg.eBOSC.pad.tfr_sample;

    for indChan = 1: numel(cfg.eBOSC.channel)
    
        display(['Channel ',num2str(indChan), '/', num2str(numel(cfg.eBOSC.channel)),': chanID ', num2str(cfg.eBOSC.channel(indChan))])
        
        cfg.tmp.channel = indChan; % encode current channel for later

        %% Step 1: time-frequency wavelet decomposition for whole signal to prepare background fit

        TFR = [];
        for indTrial = 1:numel(cfg.eBOSC.trial)
            TFR.trial{indTrial} = BOSC_tf(data.trial{indTrial}(cfg.eBOSC.channel(indChan),:),cfg.eBOSC.F,cfg.eBOSC.fsample,cfg.eBOSC.wavenumber);
        end; clear indTrial

        %% Step 2: robust background power fit (see 2020 NeuroImage paper)

        [eBOSC, pt, dt] = eBOSC_getThresholds(cfg, TFR, eBOSC);

        %% application of thresholds to single trials

        for indTrial = 1:numel(cfg.eBOSC.trial)

            cfg.tmp.trial = cfg.eBOSC.trial(indTrial); % encode current trial for later

            % get wavelet transform for single trial
            % tfr padding is removed to avoid edge artifacts from the wavelet
            % transform. Note that a padding fpr detection remains attached so that there
            % is no problems with too few sample points at the edges to
            % fulfill the duration criterion.
            TFR_ = TFR.trial{indTrial}(:,cfg.eBOSC.pad.tfr_sample+1:end-cfg.eBOSC.pad.tfr_sample);

            %% Step 3: detect rhythms and calculate Pepisode

            % The next section applies both the power and the duration
            % threshold to detect individual rhythmic segments in the continuous signals.
            detected = zeros(size(TFR_));
            for f = 1:length(cfg.eBOSC.F)
                detected(f,:) = BOSC_detect(TFR_(f,:),pt(f),dt(f),cfg.eBOSC.fsample);
            end; clear f

            % remove padding for detection (matrix with padding required for refinement)
            eBOSC.detected(indChan, indTrial,:,:) = detected(:,cfg.eBOSC.pad.detection_sample+1:end-cfg.eBOSC.pad.detection_sample);

            % encode pepisode of detected rhythms (optional)
            eBOSC.pepisode(indChan, indTrial,:) = mean(eBOSC.detected(indChan, indTrial,:,:),4);
            
            %% Step 4 (optional): create table of separate rhythmic episodes

            cfg.tmp.inputTime = data.time{cfg.tmp.trial};
            cfg.tmp.detectedTime = cfg.tmp.inputTime(cfg.eBOSC.pad.tfr_sample+1:end-cfg.eBOSC.pad.tfr_sample);
            cfg.tmp.finalTime = cfg.tmp.inputTime(cfg.eBOSC.pad.total_sample+1:end-cfg.eBOSC.pad.total_sample);
            
            [eBOSC.episodes, detected_ep] = eBOSC_episode_create(cfg,TFR_,detected,eBOSC);
            
            % remove padding for detection (already done for eBOSC.episodes)
            eBOSC.detected_ep(indChan, indTrial,:,:) = detected_ep(:,cfg.eBOSC.pad.detection_sample+1:end-cfg.eBOSC.pad.detection_sample);
            clear detected_ep;

            % encode abundance of eBOSC.episodes (optional)
            eBOSC.abundance_ep(indChan, indTrial,:) = mean(squeeze(eBOSC.detected_ep(indChan, indTrial,:,:)),2);

%           % Supplementary Plot: original eBOSC.detected vs. sparse episode power
%           figure; 
%           subplot(121); imagesc(squeeze(eBOSC.detected).*TFR_(:,cfg.eBOSC.pad.detection_sample+1:end-cfg.eBOSC.pad.detection_sample));
%           subplot(122); imagesc(squeeze(eBOSC.detected_ep).*TFR_(:,cfg.eBOSC.pad.detection_sample+1:end-cfg.eBOSC.pad.detection_sample));

        end; clear indTrial; % trial loop

    end % channel loop
end