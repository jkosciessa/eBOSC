function [eBOSC, cfg] = eBOSC_wrapper(cfg, data, eBOSC)


    % check field entries
    
    % get dimensions
    
    % initialize matrices

    
    for indChan = 1: numel(cfg.eBOSC.channel)
    
        e = cfg.eBOSC.channel(indChan);

        display(['channel #' num2str(e)])

        cfg.tmp.inputTime = data.time{1,1};
        cfg.tmp.detectedTime = cfg.tmp.inputTime(cfg.eBOSC.pad.tfr_sample+1:end-cfg.eBOSC.pad.tfr_sample);
        cfg.tmp.finalTime = cfg.tmp.inputTime(cfg.eBOSC.pad.total_sample+1:end-cfg.eBOSC.pad.total_sample);
        cfg.tmp.channel = e; % encode current channel for later

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
                % These are two alternative ways to extract the onset timepoint
                % from the table
                idx_onsetTime(indEp) = find(cfg.tmp.finalTime>= eBOSC.episodes.Onset(idx_alpha(indEp)), 1, 'first');
                idx_onset(indEp) = eBOSC.episodes.ColID{idx_alpha(indEp)}(1);
            end

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

        end; clear indTrial; % trial loop

    end % channel loop
end