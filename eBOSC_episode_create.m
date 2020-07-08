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
%    Copyright 2018 Julian Q. Kosciessa, Thomas H. Grandy, Douglas D. Garrett & Markus Werkle-Bergner.

function [detected_new,episodes] = eBOSC_episode_create(B,detected,cfg)

% This function creates continuous rhythmic "episodes" and attempts to control for the impact of wavelet parameters.
%  Time-frequency points that best represent neural rhythms are identified by
%  heuristically removing temporal and frequency leakage. 
%
% Frequency leakage: at each frequency x time point, power has to exceed neighboring frequencies.
%  Then it is checked whether the detected time-frequency points belong to
%  a continuous episode for which (1) the frequency maximally changes by 
%  +/- n steps (cfg.eBOSC.fstp) from on time point to the next and (2) that is at 
%  least as long as n number of cycles (cfg.eBOSC.ncyc) of the average freqency
%  of that episode (a priori duration threshold).
%
% Temporal leakage: The impact of the amplitude at each time point within a rhythmic episode on previous
%  and following time points is tested with the goal to exclude supra-threshold time
%  points that are due to the wavelet extension in time. 
%
%  input:   B        = time-frequency matrix
%           detected = detected oscillations in B (based on power and duration threshold)
%           cfg      - .eBOSC.F           = frequency resolution of B (log-scaled!)
%                    - .eBOSC.fsample     = sampling frequency
%                    - .eBOSC.wavenumber  = wavenumber in time-frequency analysis
%                    - .eBOSC.npnts       = length of data (in data points)
%                    - .eBOSC.fstp        = maximal step size of frequencies from one
%                                     time point to the next (steps refer
%                                     to vector entries in .F)
%                    - .eBOSC.ncyc        = minimum number of cycles of the average 
%                                     freqency of that segment
%                    - .eBOSC.BiasCorrection; apply bias correction
%                    - .eBOSC.method; method to use for bias correction; options: 'MaxBias', 'FWHM'; see Kosciessa et al. supplement for more infos
%                    - .eBOSC.edgeOnly; correct for wavelet bias only at segment edges ('yes') or also within detected segments ('no')
%                    - .eBOSC.effSignal; signal to base wavelet correaction on: 'all' - entire signal including background; 'PT' - only signal above the power threshold                
%
% No default exists for any of the parameters.
%
%  output:  detected_new = new detected matrix with frequency leakage removed
%           episodes     = {nx4} with specific episode information
%                          {n,1} = episode indices
%                          {n,2} = average frequency
%                          {n,3} = episode length (in sec)
% 

%% Accounting for the frequency spread of the wavelet

% Here, we compute the bandpass response as given by the wavelet
% formula and apply half of the BP repsonse on top of the center frequency.
% Because of log-scaling, the widths are not the same on both sides.

freqWidth = (2/cfg.eBOSC.wavenumber)*cfg.eBOSC.F;
lowFreq = cfg.eBOSC.F-(freqWidth/2);
highFreq = cfg.eBOSC.F+(freqWidth/2);

indF_low = [];
indF_high = [];
for indF = 1:numel(cfg.eBOSC.F)
	if ~isempty(find(cfg.eBOSC.F<=lowFreq(indF), 1, 'last'))
		fmat(indF,1) = find(cfg.eBOSC.F<=lowFreq(indF), 1, 'last')+1; % first freq falling into range
	else fmat(indF,1) = 1;
	end
	if ~isempty(find(cfg.eBOSC.F>=highFreq(indF), 1, 'first'))
		fmat(indF,3) = find(cfg.eBOSC.F>=highFreq(indF),1, 'first')-1; % last freq falling into range
	else fmat(indF,3) = numel(cfg.eBOSC.F);
	end
end; fmat(:,2) = (1:numel(cfg.eBOSC.F))'; clear indF;
range = diff(fmat,[],2);
range = max(range,[],1);

% initialize variables
% append search space (i.e. range at both ends. first index refers to lower range
tmp_B    = [zeros(range(1,1),size(B,2)); B.*detected; zeros(range(1,2),size(B,2))];
%tmp_B    = B.*detected;
detected = zeros(size(detected));

for f = range(1,1)+1:size(tmp_B,1)-range(1,2) % loop across frequencies. note that indexing respects the appended segments

	r = 1;
	% encode detected positions at current frequency
	tmp_det(r,:) = tmp_B(f,:) ~= 0;
	% encode detected positions within LOWER range where current freq
	% point is larger
	for r = 1:range(1,1)
		tmp_det(r,:) = tmp_B(f,:) ~= 0 & tmp_B(f,:) >= tmp_B(f-r,:);
	end
	% encode detected positions within HIGHER range where current freq
	% point is larger
	for r = 1:range(1,2)
		r2 = range(1,1)+r;
		tmp_det(r2,:) = tmp_B(f,:) ~= 0 & tmp_B(f,:) >= tmp_B(f+r,:);
	end; clear r r2;

	detected(f,:) = mean(tmp_det,1) == 1;

end; clear f;

% remove padded zeros
detected = detected(range(1,1)+1:size(tmp_B,1)-range(1,2),:);

%%  Create continuous rhythmic episodes

% add zeros
detected1        = [zeros(cfg.eBOSC.fstp,size(detected,2)); detected; zeros(cfg.eBOSC.fstp,size(detected,2))];
detected1(:,1)   = 0;
detected1(:,end) = 0;
tmp_B1        = [zeros(cfg.eBOSC.fstp,size(detected,2)); B.*detected; zeros(cfg.eBOSC.fstp,size(detected,2))];
tmp_B1(:,1)   = 0;
tmp_B1(:,end) = 0;
detected2        = zeros(size(detected));
detected_new     = zeros(size(detected));

% segment counter
j = 1;

while sum(sum(detected1)) > 0
    % sampling point counter
    k = 1;
    % find seed
    [x(k),y(k)] = find(detected1==1,1);
    % check next sampling point
    chck = 0;
    while chck == 0
        % next sampling point
        tmp = find(detected1(x(k)-cfg.eBOSC.fstp:x(k)+cfg.eBOSC.fstp,y(k)+1)==1);
        if ~isempty(tmp)
            y(k+1) = y(k) + 1;
            if numel(tmp) > 1 
                % JQK 161017: It is possible that an episode is branching 
                % two ways, hence we follow the 'strongest' branch; 
                % Note that there is no correction for 1/f here, but 
                % practically, it leads to satisfying results 
                % (i.e. following the longer episodes).
                tmp_data = tmp_B1(x(k)-cfg.eBOSC.fstp:x(k)+cfg.eBOSC.fstp,y(k)+1);
                tmp = find(tmp_data == max(tmp_data));
            end
            x(k+1) = x(k) + tmp - cfg.eBOSC.fstp - 1;
            k = k + 1;
        else
            chck = 1;
        end
    end
        
    % check for passing the duration requirement
    % get average frequency
    avg_frq = mean(cfg.eBOSC.F(x-cfg.eBOSC.fstp));
    % match to closest frequency
    [~, indF] = min(abs(cfg.eBOSC.F-avg_frq));
    % check number of data points to fulfill number of cycles criterion
    num_pnt = floor((cfg.eBOSC.fsample ./ avg_frq) .* (cfg.eBOSC.ncyc(indF))); clear indF;
    
    if length(y) >= num_pnt
        episodes{j,1} = single([x'-cfg.eBOSC.fstp,y'-cfg.eBOSC.fstp]);
        for m = 1:length(y)
            episodes{j,2}(m,1) = single(cfg.eBOSC.F(episodes{j,1}(m,1)));
            episodes{j,2}(m,2) = single(B(episodes{j,1}(m,1),episodes{j,1}(m,2)));
        end
        episodes{j,3} = single(avg_frq);
        episodes{j,4} = single(length(y) ./ cfg.eBOSC.fsample);
        for l = 1:length(y)
            detected1(x(l),y(l))                   = 0;
            detected_new(episodes{j,1}(l,1),episodes{j,1}(l,2)) = 1;
        end
        j = j + 1;
    else
        for l = 1:length(y)
            detected1(x(l),y(l))                   = 0;
        end        
    end
    % clear variables
    clear k x y chck tmp avg_frq num_pnt m
end

if sum(sum(detected_new)) == 0
    episodes = {};
end

%%  Exclude temporal amplitude "leakage" due to wavelet smearing

if strcmp(cfg.eBOSC.BiasCorrection, 'yes')
    if strcmp(cfg.eBOSC.method, 'HM') % FWHM correction
        [episodes, detected_new] = eBOSC_episode_postproc_fwhm(episodes,cfg, B);
    elseif strcmp(cfg.eBOSC.method,'MaxBias')
        [episodes, detected_new] = eBOSC_episode_postproc_maxbias(episodes,cfg, B);
    end
else
    
%% remove episodes and part of episodes that fall into 'shoulder'

[episodes] = eBOSC_episode_rm_shoulder(cfg,detected_new,episodes);

end
