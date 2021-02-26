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

function [episodes_new, detected_new] = eBOSC_episode_postproc_maxbias(cfg, episodes, TFR)
% This function performs post-processing of input episodes by checking
% whether 'detected' time points can be explained by the simulated extension of
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

disp("Applying maxbias post-processing ...")

% re-initialize detected_new (for post-proc results)
N_freq = size(TFR,1); N_tp = size(TFR,2);
detected_new = zeros(N_freq, N_tp);

% check if there are episodes
if ~isempty(episodes)

    % generate "bias" matrix
    B_bias = zeros(length(cfg.eBOSC.F),numel(cfg.eBOSC.F),2*N_tp+1);

    amp_max = NaN(length(cfg.eBOSC.F), length(cfg.eBOSC.F));
    for f = 1:length(cfg.eBOSC.F)
        % temporary time vector and signal
        time = 1/cfg.eBOSC.fsample:1/cfg.eBOSC.fsample:(1/cfg.eBOSC.F(f));
        tmp_sig = cos(time*2*pi*cfg.eBOSC.F(f)).*-1+1;
        % signal for time-frequency analysis
        signal = [zeros(1,N_tp) tmp_sig zeros(1,N_tp)];
        tmp_bias_mat = BOSC_tf(signal,cfg.eBOSC.F,cfg.eBOSC.fsample,cfg.eBOSC.wavenumber);
        % bias matrix
        B_bias(f,:,1:N_tp+1)   = tmp_bias_mat(:,1:N_tp+1);
        B_bias(f,:,N_tp+1:end) = fliplr(tmp_bias_mat(:,1:N_tp+1));
        % maximum amplitude
        amp_max(f,:) = max(squeeze(B_bias(f,:,:)),[],2);
        % clear variables
        clear B_tmp ind_* signal time tmp_sig
    end; clear f

    % midpoint index
    ind_mid = N_tp+1;

    % set counter
    cnt = 1;

    % loop episodes
    for e = 1:size(episodes,1)

        % get temporary frequency vector
        f_ = episodes.Frequency{e}';
        % get temporary amplitude vector
        a_ = episodes.Power{e};

        m_ = zeros(length(a_),length(a_));
        
        % location in time with reference to matrix TFR
        t_ind = episodes.ColID{e}(1):episodes.ColID{e}(end);

        % indices of time points' frequencies within "bias" matrix
        f_vec = episodes.RowID{e};

        %figure; hold on;
        for tp = 1:length(a_)
            % index of current point's frequency within "bias" matrix
            ind_f = f_vec(tp);
            % get bias vector that varies with frequency of the
            % timepoints in the episode
            temporalBiasIndices = ind_mid+1-tp:ind_mid+length(a_)-tp;
            tmp_biasVec = B_bias(sub2ind(size(B_bias),repmat(ind_f,numel(f_vec),1),f_vec,temporalBiasIndices'))';
            % temporary "bias" vector (frequency-varying)
            if strcmp(cfg.eBOSC.postproc.effSignal,'all')
                tmp_bias = ((tmp_biasVec./amp_max(ind_f,f_vec))*a_(tp));
            elseif strcmp(cfg.eBOSC.postproc.effSignal,'PT')
                tmp_bias = ((tmp_biasVec./amp_max(ind_f,f_vec))*(a_(tp)-cfg.tmp.pt(ind_f))) + cfg.tmp.pt(ind_f);
            end
            % compare to data
            m_(tp,:) = a_' >= tmp_bias;
            %plot(a_', 'k'); hold on; plot(tmp_bias, 'r');
            % clear variables
            clear ind_f amp_fctr tmp_bias tmp_biasVec
        end; clear tp

        % identify which time points to retain and discard
        % Options: only correct at signal edge; correct within entire signal
        if strcmp(cfg.eBOSC.postproc.edgeOnly,'no')
            keep = sum(m_) == length(a_);
        elseif strcmp(cfg.eBOSC.postproc.edgeOnly,'yes')
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