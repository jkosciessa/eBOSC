function [episodes_new, detected_new] = eBOSC_episode_postproc_maxbias(episodes,cfg, B)

% This method works as follows: we estimate the bias introduced by
% wavelet convolution. The bias is represented by the amplitudes
% estimated for the zero-shouldered signal (i.e. for which no real 
% data was initially available). The influence of episodic
% amplitudes on neighboring time points is assessed by scaling each
% time point's amplitude with the last 'rhythmic simulated time
% point', i.e. the first time wavelet amplitude in the simulated
% rhythmic time points. At this time point the 'bias' is maximal,
% although more precisely, this amplitude does not represent a
% bias per se.

% check if there are episodes
if ~isempty(episodes)

    % generate "bias" matrix
    B_bias = zeros(length(cfg.eBOSC.F),numel(cfg.eBOSC.F),2*cfg.eBOSC.npnts+1);

    for f = 1:length(cfg.eBOSC.F)
        % temporary time vector and signal
        time = 1/cfg.eBOSC.fsample:1/cfg.eBOSC.fsample:(1/cfg.eBOSC.F(f));
        tmp_sig = cos(time*2*pi*cfg.eBOSC.F(f)).*-1+1;
        % signal for time-frequency analysis
        signal = [zeros(1,cfg.eBOSC.npnts) tmp_sig zeros(1,cfg.eBOSC.npnts)];
        tmp_bias_mat = BOSC_tf(signal,cfg.eBOSC.F,cfg.eBOSC.fsample,cfg.eBOSC.wavenumber);
        % bias matrix
        B_bias(f,:,1:cfg.eBOSC.npnts+1)   = tmp_bias_mat(:,1:cfg.eBOSC.npnts+1);
        B_bias(f,:,cfg.eBOSC.npnts+1:end) = fliplr(tmp_bias_mat(:,1:cfg.eBOSC.npnts+1));
        % maximum amplitude
        amp_max(f,:) = max(squeeze(B_bias(f,:,:)),[],2);
        % clear variables
        clear B_tmp ind_* signal time tmp_sig
    end; clear f

    % midpoint index
    ind_mid = cfg.eBOSC.npnts+1;

    % set counter
    cnt = 1;

    % loop episodes
    for e = 1:size(episodes,1)

        % get temporary frequency vector
        f_ = episodes{e,2}(:,1);

        % get temporary amplitude vector
        a_ = episodes{e,2}(:,2);

        m_ = zeros(length(a_),length(a_));

        % indices of time points' frequencies within "bias" matrix
        f_vec = episodes{e,1}(:,1);

        %figure; hold on;
        for tp = 1:length(episodes{e,2})
            % index of current point's frequency within "bias" matrix
            ind_f = f_vec(tp);
            % get bias vector that varies with frequency of the
            % timepoints in the episode
            temporalBiasIndices = ind_mid+1-tp:ind_mid+length(a_)-tp;
            tmp_biasVec = B_bias(sub2ind(size(B_bias),repmat(ind_f,numel(f_vec),1),f_vec,temporalBiasIndices'))';
            % temporary "bias" vector (frequency-varying)
            if strcmp(cfg.eBOSC.effSignal,'all')
                tmp_bias = ((tmp_biasVec./amp_max(ind_f,f_vec))*a_(tp));
            elseif strcmp(cfg.eBOSC.effSignal,'PT')
                tmp_bias = ((tmp_biasVec./amp_max(ind_f,f_vec))*(a_(tp)-cfg.eBOSC.pt(ind_f))) + cfg.eBOSC.pt(ind_f);
            end
            % compare to data
            m_(tp,:) = a_' >= tmp_bias;
            %plot(a_', 'k'); hold on; plot(tmp_bias, 'r');
            % clear variables
            clear ind_f amp_fctr tmp_bias tmp_biasVec
        end; clear tp

        % identify which time points to retain and discard
        % Options: only correct at signal edge; correct within entire signal
        if strcmp(cfg.eBOSC.edgeOnly,'no')
            keep = sum(m_) == length(a_);
        elseif strcmp(cfg.eBOSC.edgeOnly,'yes')
            keep = sum(m_) == length(a_);
            keepEdgeRemovalOnly = zeros(1, numel(keep));
            keepEdgeRemovalOnly(find(keep==1,1,'first'):find(keep==1,1,'last')) = 1;
            keep = keepEdgeRemovalOnly;
        end

        % get new episodes
        keep   = [0 keep 0];
        d_keep = diff(keep);

        if max(d_keep) == 1 && min(d_keep) == -1

        % start and end indices
        ind_epsd(:,1) = find(d_keep ==  1);
        ind_epsd(:,2) = find(d_keep == -1)-1;

        for i = 1:size(ind_epsd,1)
            % temporary frequency & amplitude vector
            tmp = episodes{e,2}(ind_epsd(i,1):ind_epsd(i,2),:);
            % check for passing the duration requirement
            % get average frequency
            avg_frq = mean(tmp(:,1));
            % match to closest frequency
            [~, indF] = min(abs(cfg.eBOSC.F-avg_frq));
            % check number of data points to fulfill number of cycles criterion
            num_pnt = floor((cfg.eBOSC.fsample ./ avg_frq) .* (cfg.eBOSC.ncyc(indF))); clear indF;
            if num_pnt <= size(tmp,1)
                episodes_new{cnt,1} = episodes{e,1}(ind_epsd(i,1):ind_epsd(i,2),:);
                episodes_new{cnt,2} = episodes{e,2}(ind_epsd(i,1):ind_epsd(i,2),:);
                episodes_new{cnt,3} = single(avg_frq);
                episodes_new{cnt,4} = single(size(tmp,1) ./ cfg.eBOSC.fsample);        
                for l = 1:size(tmp,1)
                    detected_new(episodes_new{cnt,1}(l,1),episodes_new{cnt,1}(l,2)) = 1;
                end; clear l
                % set counter
                cnt = cnt + 1;
            end
            % clear variables
            clear avg_frq num_pnt tmp
        end
        % clear variables
        clear ind_epsd
        end
        % clear variables
        clear f_ a_ m_ keep d_keep
    end; clear e
    
end