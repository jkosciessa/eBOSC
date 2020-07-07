% check stabilty of single-trial abundance, amplitude

% 170830: initially wrong referencing of sessions (e.g. CuttingEEG); now corrected
% 181002: average across loads

% pn.root     = '/Volumes/EEG/BOSC_Sternberg/B_analyses/B_BOSC/B_eBOSC_170703_full_wl6/';
% pn.tools    = [pn.root, 'C_interrelations/T_tools/']; addpath(genpath(pn.tools));
% pn.data     = [pn.root, 'B_extractIndices/B_data/'];
% 
% load([pn.data, 'X_8to15_170703_v3.mat']);
% 
% %% create data for plot
% 
% Figure6C.info.channels = 44:60;
% Figure6C.info.cond = {'L2','L4','L6'};
% Figure6C.info.sessions = {'1', '7', '8'};
% Figure6C.info.colorm = brewermap(4,'greys'); % use external function to create colormap
% 
% for indSession = 1:3
%     abundance = cat(4, X{indSession,1}.L2e_abn, X{indSession,1}.L4e_abn, X{indSession,1}.L6e_abn); % sub*chan*trial
%     abundance = squeeze(nanmean(abundance,4));
%     chanAvg = squeeze(nanmean(abundance(:,Figure6C.info.channels,:),2));
%     % sort with reference to S1
%     if indSession == 1
%         chanAvgToSort = chanAvg;
%     end
%     [~, index] = sort(nanmean(chanAvgToSort,2));
%     % get average
%     meanSorted = nanmean(chanAvg,2);
%     Figure6C.dataToPlot.meanSorted(indSession,:) = meanSorted(index);
%     % get STD
%     stdSorted = nanstd(chanAvg,[],2);
%     Figure6C.dataToPlot.stdSorted(indSession,:) = stdSorted(index);
% end
% 
% %% save Figure data
% 
% save('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/F6C.mat', 'Figure6C')

%% load Figure data

load('/Users/kosciessa/Desktop/eBOSC/figureData/F5C.mat', 'Figure6C')

%% FIGURE 3A: Individual abundance stability

h = figure('units','normalized','position',[.1 .1 .7 .2]);
for indSession = 1:3
    if indSession == 1
        x = .75:1:31.75;
        design = 'sk';
    elseif indSession == 2
        x = 1:1:32;
        design = 'sr';
    elseif indSession == 3
        x = 1.25:1:32.25;
        design = 'sb';
    end
    e = errorbar(x, Figure6C.dataToPlot.meanSorted(indSession,:), Figure6C.dataToPlot.stdSorted(indSession,:), design); hold on;
    e.Marker = 'o';
    e.MarkerSize = 10;
    e.Color = Figure6C.info.colorm(indSession+1,:);
    e.MarkerEdgeColor = Figure6C.info.colorm(indSession+1,:);
    e.MarkerFaceColor = Figure6C.info.colorm(indSession+1,:);
end
title(['Individual abundance stability']);
xlabel('Subject (sorted by mean abundance)'); ylabel({'Abundance'; '(Mean +-STD)'});
legend({'Session 1'; 'Session 7'; 'Session 8'}, 'Orientation', 'horizontal', 'location', 'NorthWest');
legend('boxoff');
xlim([.5 32.5]); ylim([-.1 1]);
xticks([1:32])
xticklabels({1:32})
ax = gca;
ax.TickDir = 'out';
set(gca,'box','off')

set(findall(gcf,'-property','FontSize'),'FontSize',18)

%% correlation matrix of Figure6C.info.sessions [reported in manuscript]

[rMat, pMat] = corrcoef(Figure6C.dataToPlot.meanSorted')

% rMat =
% 
%     1.0000    0.9494    0.9299
%     0.9494    1.0000    0.9608
%     0.9299    0.9608    1.0000

