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

function [episodes] = eBOSC_episode_rm_shoulder(cfg,detected1,episodes)
% Remove parts of the episode that fall into the 'shoulder' of individual
% trials.
% There is no check for adherence to a given duration criterion necessary,
% as the point of the padding of the detected matrix is exactly to account
% for allowing the presence of a few cycles.

    ind1 = cfg.eBOSC.pad.detection_sample+1;
    ind2 = size(detected1,2) - cfg.eBOSC.pad.detection_sample;
    cnt = 1;
    rmv = [];
    for j = 1:size(episodes,1)
        % get time points of current episode
        tmp_col = episodes.ColID{j}(1):episodes.ColID{j}(2);
        % find time points that fall inside the padding (i.e. on- and offset)
        ex = find(tmp_col < ind1 | tmp_col > ind2);
        % remove padded time points from episodes
        tmp_col(ex) = [];
        episodes.RowID{j}(ex) = [];
        episodes.Power{j}(ex) = [];
        episodes.Frequency{j}(ex) = [];
        episodes.SNR{j}(ex) = [];
        % if nothing remains of episode: track for later deletion
        if isempty(tmp_col)
            rmv(cnt,1) = j;
            cnt = cnt + 1;
            % skip further encoding to avoid problems
            continue;
        end
        % shift onset according to padding
        % Important: new col index is indexing w.r.t. to matrix AFTER
        % detected padding is removed!
        tmp_col = tmp_col - ind1 + 1;
        episodes.ColID{j} = [tmp_col(1), tmp_col(end)]';
        % WIP: episodes{j,1}(:,2) = episodes.ColID{j}
        % re-compute mean frequency
        episodes.FrequencyMean(j) = mean(episodes.Frequency{j});
        % re-compute mean amplitude
        episodes.PowerMean(j) = mean(episodes.Power{j});
        % re-compute mean SNR
        episodes.SNRMean(j) = nanmean(episodes.SNR{j});
        % re-compute duration
        episodes.DurationS(j) = (episodes.ColID{j}(2)-episodes.ColID{j}(1)+1) / cfg.eBOSC.fsample; clear ex
        episodes.DurationC(j) = episodes.DurationS(j)*episodes.FrequencyMean(j);
        % update absolute on-/offsets
        episodes.Onset(j) = cfg.tmp.finalTime(episodes.ColID{j}(1)); % episode onset in absolute time
        episodes.Offset(j) = cfg.tmp.finalTime(episodes.ColID{j}(end)); % episode offset in absolute time
    end; clear j cnt
    % remove now empty episodes from table    
    episodes(rmv,:) = []; clear rmv

end % function end