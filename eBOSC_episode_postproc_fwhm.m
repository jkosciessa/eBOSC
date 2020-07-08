function [episodes_new] = eBOSC_episode_postproc_fwhm(episodes,cfg, B)

if ~isempty(episodes)
            
    % set counter
    cnt = 1;

    for e = 1:size(episodes,1)

        % get temporary frequency vector
        f_ = episodes{e,2}(:,1);
        f_unique = unique(f_);
        f_ind_unique = find(ismember(cfg.eBOSC.F', f_unique, 'rows'));

        % get temporary amplitude vector
        a_ = episodes{e,2}(:,2);

        % location in time with reference to matrix B
        t_ind = episodes{e,1}(:,2);

        % initiate bias matrix
        biasMat = [];
        biasMat = zeros(numel(f_unique),numel(a_));

        for tp = 1:numel(a_)

            % The FWHM correction is done independently at each
            % frequency. To accomplish this, we actually reference
            % to the original data in the TF matrix.

            for f = 1:numel(f_unique) % search within frequencies that occur within the episode
                % create wavelet with center frequency and amplitude at time point
                st=1./(2*pi*(f_unique/cfg.eBOSC.wavenumber));
                t=-3.6*st(f):(1/cfg.eBOSC.fsample):3.6*st(f);
                if strcmp(cfg.eBOSC.effSignal, 'all')
                    m{f}=(B(f_ind_unique(f), t_ind(tp)))*exp(-t.^2/(2*st(f)^2)).*exp(1i*2*pi*f_unique(f).*t); % Morlet wavelet with amplitude-power threshold modulation
                elseif strcmp(cfg.eBOSC.effSignal, 'PT')
                    m{f}=(B(f_ind_unique(f), t_ind(tp))-cfg.eBOSC.pt(f_ind_unique(f)))*exp(-t.^2/(2*st(f)^2)).*exp(1i*2*pi*f_unique(f).*t); % Morlet wavelet with amplitude-power threshold modulation
                end
                wl_a = []; wl_a = abs(m{f}); % amplitude of wavelet
                [maxval(f), maxloc(f)] = max(wl_a);
                index_fwhm(f) = find(wl_a>= maxval(f)/2, 1, 'first');
                fwhm_a(f) = wl_a(index_fwhm(f)); % amplitude at fwhm, freq
                if strcmp(cfg.eBOSC.effSignal, 'PT')
                    fwhm_a(f) = fwhm_a(f)+cfg.eBOSC.pt(f_ind_unique(f)); % re-add power threshold
                end
                correctionDist(f) = maxloc(f)-index_fwhm(f);
                % extract FWHM amplitude of frequency- and amplitude-specific wavelet
                if tp-correctionDist(f) > 0 % check that lower fwhm is part of signal 
                    if biasMat(f,tp-correctionDist(f)) < fwhm_a(f) % and that existing value is lower than update
                        biasMat(f,tp-correctionDist(f)) = fwhm_a(f);
                    end
                end
                if tp+correctionDist(f) <= size(biasMat,2) % check that upper fwhm is part of signal 
                    if biasMat(f,tp+correctionDist(f)) < fwhm_a(f) % and that existing value is lower than update
                        biasMat(f,tp+correctionDist(f)) = fwhm_a(f);
                    end
                end 
            end % frequencies within episode
        end % timepoints

        % retain only those points that are larger than the FWHM
        aMat_retain = [];
        [~, indFreqs] = ismember(f_, f_unique);
        for indF = 1:numel(f_unique)
             aMat_retain(indF,indFreqs == indF) = a_(indFreqs == indF)';
        end
        aMat_retain(aMat_retain <= biasMat) = 0; % anything that is lower than the convolved wavelet is removed

        % identify which time points to retain and discard
        % Options: only correct at signal edge; correct within entire signal
        if strcmp(cfg.eBOSC.edgeOnly,'no')
            keep = mean(aMat_retain,1)>0;
            keep = (keep)>0;
        elseif strcmp(cfg.eBOSC.edgeOnly,'yes')
            keep = mean(aMat_retain,1)>0;
            keep = (keep)>0;
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