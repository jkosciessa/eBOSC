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

    ind1 = cfg.eBOSC.detectedPad+1;
    ind2 = size(detected1,2) - cfg.eBOSC.detectedPad;
    cnt = 1;
    rmv = [];
    for j = 1:size(episodes,1)
        % final episodes
        ex = find(episodes{j,1}(:,2) < ind1 | episodes{j,1}(:,2) > ind2);
        % update episodes
        episodes{j,1}(ex,:) = [];
        episodes{j,1}(:,2) = episodes{j,1}(:,2) - ind1 + 1;
        episodes{j,2}(ex,:) = [];
        episodes{j,3}       = mean(episodes{j,2}(:,1));
        episodes{j,4}       = size(episodes{j,1},1) / cfg.eBOSC.fsample; clear ex
        if isempty(episodes{j,1})
            rmv(cnt,1) = j;
            cnt = cnt + 1;
        end
        % original episodes
        ex = find(episodes{j,5}(:,2) < ind1 | episodes{j,5}(:,2) > ind2);
        % update episodes
        episodes{j,5}(ex,:) = [];
        episodes{j,5}(:,2) = episodes{j,5}(:,2) - ind1 + 1;
        episodes{j,6}(ex,:) = [];
        episodes{j,7}       = mean(episodes{j,6}(:,1));
        episodes{j,8}       = size(episodes{j,5},1) / cfg.eBOSC.fsample; clear ex
    end; clear j cnt
    episodes(rmv,:) = []; clear rmv


end % function end