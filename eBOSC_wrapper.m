function [eBOSC, cfg] = eBOSC_wrapper(cfg, data)

    eBOSC = [];
    
    % set some defaults for included channels and trials, if not specified
    if isempty(cfg.eBOSC.channel)
        cfg.eBOSC.channel = 1:numel(data.label);
    end
    if isempty(cfg.eBOSC.trial)
        cfg.eBOSC.trial = 1:numel(data.trial);
    end
    if isempty(cfg.eBOSC.trial)
        cfg.eBOSC.trial_background = 1:numel(data.trial);
    end

    for indChan = 1: numel(cfg.eBOSC.channel)
    
        e = cfg.eBOSC.channel(indChan);

        display(['Channel ',num2str(indChan), '/', num2str(numel(cfg.eBOSC.channel)),': chanID ', num2str(e)])

        cfg.tmp.inputTime = data.time{1,1};
        cfg.tmp.detectedTime = cfg.tmp.inputTime(cfg.eBOSC.pad.tfr_sample+1:end-cfg.eBOSC.pad.tfr_sample);
        cfg.tmp.finalTime = cfg.tmp.inputTime(cfg.eBOSC.pad.total_sample+1:end-cfg.eBOSC.pad.total_sample);
        cfg.tmp.channel = [indChan, e]; % encode current channel for later

        %% Step 1: time-frequency wavelet decomposition for whole signal to prepare background fit

        eBOSC.Ntrial = length(cfg.eBOSC.trial);

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

        %% application of thresholds to single trials

        for indTrial = 1:numel(cfg.eBOSC.trial)

            cfg.tmp.trial = cfg.eBOSC.trial(indTrial); % encode current trial for later

            % get wavelet transform for single trial
            % tfr padding is removed to avoid edge artifacts from the wavelet
            % transform. Note that a padding fpr detection remains attached so that there
            % is no problems with too few sample points at the edges to
            % fulfill the numcycles criterion.
            TFR_ = TFR.trial{1,indTrial}(:,cfg.eBOSC.pad.tfr_sample+1:end-cfg.eBOSC.pad.tfr_sample);

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

            % encode original signals (optional)
            origData = data.trial{indTrial}(e, cfg.eBOSC.pad.total_sample+1:end-cfg.eBOSC.pad.total_sample);
            eBOSC.origData(indChan, indTrial,:) = origData;
            
            %% Step 4 (optional): create table of separate rhythmic episodes

            [detected_ep,eBOSC.episodes] = eBOSC_episode_create(TFR_,eBOSC, cfg, detected);

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