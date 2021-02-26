%    This file is part of the extended Better OSCillation detection (eBOSC) library.
%
%    The eBOSC library is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    The eBOSC library is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
%
%    Copyright 2020 Julian Q. Kosciessa, Thomas H. Grandy, Douglas D. Garrett & Markus Werkle-Bergner

function [episodes_new, detected_new] = eBOSC_episode_postproc_fwhm(cfg, episodes, TFR)
% This function performs post-processing of input episodes by checking
% whether 'detected' time points can trivially be explained by the FWHM of
% the wavelet used in the time-frequency transform.
%
% Inputs: 
%           cfg | config structure with cfg.eBOSC field
%           episodes | table of episodes
%           TFR | time-frequency matrix
%
% Outputs: 
%           episodes_new | updated table of episodes
%           detected_new | updated binary detected matrix

disp("Applying FWHM post-processing ...")
% re-initialize detected_new (for post-proc results)
detected_new = zeros(size(TFR));

if ~isempty(episodes)
            
    % set counter
    cnt = 1;

    for e = 1:size(episodes,1)
        
        % get temporary frequency vector
        f_ = episodes.Frequency{e}';
        f_unique = unique(f_);
        f_ind_unique = find(ismember(cfg.eBOSC.F', f_unique', 'rows'));

        % get temporary amplitude vector
        a_ = episodes.Power{e};

        % location in time with reference to matrix TFR
        t_ind = episodes.ColID{e}(1):episodes.ColID{e}(end);

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
                if strcmp(cfg.eBOSC.postproc.effSignal, 'all')
                    m{f}=(TFR(f_ind_unique(f), t_ind(tp)))*exp(-t.^2/(2*st(f)^2)).*exp(1i*2*pi*f_unique(f).*t); % Morlet wavelet with amplitude-power threshold modulation
                elseif strcmp(cfg.eBOSC.postproc.effSignal, 'PT')
                    m{f}=(TFR(f_ind_unique(f), t_ind(tp))-cfg.tmp.pt(f_ind_unique(f)))*exp(-t.^2/(2*st(f)^2)).*exp(1i*2*pi*f_unique(f).*t); % Morlet wavelet with amplitude-power threshold modulation
                end
                wl_a = []; wl_a = abs(m{f}); % amplitude of wavelet
                [maxval(f), maxloc(f)] = max(wl_a);
                index_fwhm(f) = find(wl_a>= maxval(f)/2, 1, 'first');
                fwhm_a(f) = wl_a(index_fwhm(f)); % amplitude at fwhm, freq
                if strcmp(cfg.eBOSC.postproc.effSignal, 'PT')
                    fwhm_a(f) = fwhm_a(f)+cfg.tmp.pt(f_ind_unique(f)); % re-add power threshold
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
        if strcmp(cfg.eBOSC.postproc.edgeOnly,'no')
            keep = mean(aMat_retain,1)>0;
            keep = (keep)>0;
        elseif strcmp(cfg.eBOSC.postproc.edgeOnly,'yes')
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
            % check for passing the duration requirement
            % get average frequency
            tmp_col = ind_epsd(i,1):ind_epsd(i,2);
            avg_frq = mean(f_(tmp_col));
            % match to closest frequency
            [~, indF] = min(abs(cfg.eBOSC.F-avg_frq));
            % check number of data points to fulfill number of cycles criterion
            num_pnt = floor((cfg.eBOSC.fsample ./ avg_frq) .* (cfg.eBOSC.threshold.duration(indF))); clear indF;
            % if this duration criterion is still fulfilled, encode in table
            if num_pnt <= size(tmp_col,2)
                % exchange x and y with relevant info
                % update all data in table with new episode limits
                epData.row(cnt) = {episodes.RowID{e}(ind_epsd(i,1):ind_epsd(i,2))};
                epData.col(cnt) = {[t_ind(ind_epsd(i,1)), t_ind(ind_epsd(i,2))]'};
                epData.freq(cnt) = {f_(ind_epsd(i,1):ind_epsd(i,2))'};
                epData.freqMean(cnt) = single(avg_frq);
                epData.pow(cnt) = {a_(ind_epsd(i,1):ind_epsd(i,2))};
                epData.powMean(cnt) = nanmean(epData.pow{cnt});
                epData.durS(cnt) = single(length(epData.pow{cnt}) ./ cfg.eBOSC.fsample);
                epData.durC(cnt) = epData.durS(cnt)*epData.freqMean(cnt);
                epData.trial(cnt) = cfg.tmp.trial;
                epData.chan(cnt) = cfg.eBOSC.channel(cfg.tmp.channel);
                epData.onset(cnt) = cfg.tmp.detectedTime(epData.col{cnt}(1)); % episode onset in absolute time
                epData.offset(cnt) = cfg.tmp.detectedTime(epData.col{cnt}(end)); % episode offset in absolute time
                epData.snr(cnt) = {episodes.SNR{e}(ind_epsd(i,1):ind_epsd(i,2))};
                epData.snrMean(cnt) = nanmean(epData.snr{cnt});
                % set all detected points to one in binary detected matrix
                detected_new(sub2ind(size(TFR),epData.row{cnt},[epData.col{cnt}(1):epData.col{cnt}(2)]')) = 1;
                % set counter
                cnt = cnt + 1;
            end
            % clear variables
            clear avg_frq num_pnt
        end
        % clear variables
        clear ind_epsd
        end
        % clear variables
        clear f_ a_ m_ keep d_keep
    end; clear e
    
    % prepare for the contingency that no episodes are created
    varNames = {'Trial', 'Channel', 'FrequencyMean', 'DurationS', 'DurationC', 'PowerMean', 'Onset', 'Offset', 'Power', 'Frequency', 'RowID', 'ColID', 'SNR', 'SNRMean'};
    if exist('epData', 'var')
        episodes_new = table(epData.trial', epData.chan', epData.freqMean', epData.durS',epData.durC',  epData.powMean', epData.onset', epData.offset', epData.pow', epData.freq', epData.row', epData.col', epData.snr', epData.snrMean',  ...
                'VariableNames', varNames);
    else
        episodes_new  = cell2table(cell(0,numel(varNames)), 'VariableNames', varNames);
    end; clear varNames
    
end