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

function [episodesTable, detected_new] = eBOSC_episode_create(cfg,TFR,detected,eBOSC)

% This function creates continuous rhythmic "episodes" and attempts to control for the impact of wavelet parameters.
%  Time-frequency points that best represent neural rhythms are identified by
%  heuristically removing temporal and frequency leakage. 
%
% Frequency leakage: at each frequency x time point, power has to exceed neighboring frequencies.
%  Then it is checked whether the detected time-frequency points belong to
%  a continuous episode for which (1) the frequency maximally changes by 
%  +/- n steps (cfg.eBOSC.fstp) from on time point to the next and (2) that is at 
%  least as long as n number of cycles (cfg.eBOSC.threshold.duration) of the average freqency
%  of that episode (a priori duration threshold).
%
% Temporal leakage: The impact of the amplitude at each time point within a rhythmic episode on previous
%  and following time points is tested with the goal to exclude supra-threshold time
%  points that are due to the wavelet extension in time. 
%
%  Input:   
%           cfg         | config structure with cfg.eBOSC field
%           TFR         | time-frequency matrix (excl. WLpadding)
%           detected    | detected oscillations in TFR (based on power and duration threshold)
%           eBOSC       | main eBOSC output structure; necessary to read in
%                           prior eBOSC.episodes if they exist in a loop

%  Output:  
%           detected_new    | new detected matrix with frequency leakage removed
%           episodesTable   | table with specific episode information:
%                 Trial: trial index (corresponds to cfg.eBOSC.trial)
%                 Channel: channel index
%                 FrequencyMean: mean frequency of episode (Hz)
%                 DurationS: episode duration (in sec)
%                 DurationC: episode duration (in cycles, based on mean frequency)
%                 PowerMean: mean amplitude of amplitude
%                 Onset: episode onset in s
%                 Offset: episode onset in s
%                 Power: (cell) time-resolved wavelet-based amplitude estimates during episode
%                 Frequency: (cell) time-resolved wavelet-based frequency
%                 RowID: (cell) row index (frequency dimension): following eBOSC_episode_rm_shoulder relative to data excl. detection padding
%                 ColID: (cell) column index (time dimension)
%                 SNR: (cell) time-resolved signal-to-noise ratio: momentary amplitude/static background estimate at channel*frequency
%                 SNRMean: mean signal-to-noise ratio

%% Accounting for the frequency spread of the wavelet

% Here, we compute the bandpass response as given by the wavelet
% formula and apply half of the BP repsonse on top of the center frequency.
% Because of log-scaling, the widths are not the same on both sides.

detected = eBOSC_episode_sparsefreq(cfg, detected, TFR);

%%  Create continuous rhythmic episodes

cfg.eBOSC.fstp = 1; % define step size in adjacency matrix

% add zeros
detected_remaining        = [zeros(cfg.eBOSC.fstp,size(detected,2)); detected; zeros(cfg.eBOSC.fstp,size(detected,2))];
detected_remaining(:,1)   = 0;
detected_remaining(:,end) = 0;
% detected_remaining serves as a dummy matrix; unless all entries from detected_remaining are
% removed, we will continue extracting episodes
tmp_B1        = [zeros(cfg.eBOSC.fstp,size(detected,2)); TFR.*detected; zeros(cfg.eBOSC.fstp,size(detected,2))];
tmp_B1(:,1)   = 0;
tmp_B1(:,end) = 0;
detected_new     = zeros(size(detected));

% segment counter
j = 1;

while sum(sum(detected_remaining)) > 0
    % sampling point counter
    k = 1;
    % find seed
    [x(k),y(k)] = find(detected_remaining==1,1);
    % check next sampling point
    chck = 0;
    while chck == 0
        % next sampling point
        tmp = find(detected_remaining(x(k)-cfg.eBOSC.fstp:x(k)+cfg.eBOSC.fstp,y(k)+1)==1);
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
                % if there is more than one maximum, choose first solution
                if size(tmp,1)>1
                    tmp=tmp(1,1);
                end
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
    num_pnt = floor((cfg.eBOSC.fsample ./ avg_frq) .* (cfg.eBOSC.threshold.duration(indF))); clear indF;
    
    if length(y) >= num_pnt

        epData.row(j) = {single(x'-cfg.eBOSC.fstp)};
        epData.col(j) = {[single(y(1)-cfg.eBOSC.fstp), single(y(end)-cfg.eBOSC.fstp)]'};
        epData.freq(j) = {single(cfg.eBOSC.F(epData.row{j}))'};
        epData.freqMean(j) = single(avg_frq);
        epData.pow(j) = {single(TFR(sub2ind(size(TFR),epData.row{j},[epData.col{j}(1):epData.col{j}(2)]')))};
        epData.powMean(j) = nanmean(epData.pow{j});
        epData.durS(j) = single(length(y) ./ cfg.eBOSC.fsample);
        epData.durC(j) = epData.durS(j)*epData.freqMean(j);
        epData.trial(j) = cfg.tmp.trial;
        epData.chan(j) = cfg.eBOSC.channel(cfg.tmp.channel);
        epData.onset(j) = cfg.tmp.detectedTime(epData.col{j}(1)); % episode onset in absolute time
        epData.offset(j) = cfg.tmp.detectedTime(epData.col{j}(end)); % episode offset in absolute time
        epData.snr(j) = {epData.pow{j}./eBOSC.static.pt(cfg.tmp.channel(1),epData.row{j})'}; % extract (static) background power at frequencies
        epData.snrMean(j) = nanmean(epData.snr{j});
        
        for l = 1:length(y)
            detected_remaining(x(l),y(l)) = 0;
        end
        % set all detected points to one in binary detected matrix
        detected_new(sub2ind(size(TFR),epData.row{j},[epData.col{j}(1):epData.col{j}(2)]')) = 1;
        
        j = j + 1;
    else
        for l = 1:length(y)
            detected_remaining(x(l),y(l)) = 0;
        end        
    end
    % clear variables
    clear k x y chck tmp avg_frq num_pnt m
end; clear detected_remaining

% prepare for the contingency that no episodes are created
varNames = {'Trial', 'Channel', 'FrequencyMean', 'DurationS', 'DurationC', 'PowerMean', 'Onset', 'Offset', 'Power', 'Frequency', 'RowID', 'ColID', 'SNR', 'SNRMean'};
if exist('epData', 'var')
    episodesTable = table(epData.trial', epData.chan', epData.freqMean', epData.durS',epData.durC',  epData.powMean', epData.onset', epData.offset', epData.pow', epData.freq', epData.row', epData.col', epData.snr', epData.snrMean',  ...
            'VariableNames', varNames);
else
    episodesTable  = cell2table(cell(0,numel(varNames)), 'VariableNames', varNames);
end; clear varNames

%%  Exclude temporal amplitude "leakage" due to wavelet smearing

cfg.tmp.pt = eBOSC.static.pt(cfg.tmp.channel,:); % temporarily pass on power threshold 

if strcmp(cfg.eBOSC.postproc.use, 'yes') && exist('epData', 'var') % only do this if there are any episodes to fine-tune
    if strcmp(cfg.eBOSC.postproc.method, 'FWHM')
        [episodesTable, detected_new] = eBOSC_episode_postproc_fwhm(cfg, episodesTable, TFR);
    elseif strcmp(cfg.eBOSC.postproc.method,'MaxBias')
        [episodesTable, detected_new] = eBOSC_episode_postproc_maxbias(cfg, episodesTable, TFR);
    end
end
    
%% remove episodes and part of episodes that fall into 'shoulder'

if ~isempty(episodesTable)
    [episodesTable] = eBOSC_episode_rm_shoulder(cfg,detected_new,episodesTable);
end

%% if an episode list already exists, append results

if isfield(eBOSC, 'episodes')
    episodesTable = cat(1,eBOSC.episodes, episodesTable);
end

end
