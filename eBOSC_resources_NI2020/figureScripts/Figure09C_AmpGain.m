% - used to be /Volumes/EEG/BOSC_SternRest/B_scripts_JQK/B_processing_170704/E_JQKplotAmpAbn_170816B.mlx

% pn.root = '/Volumes/EEG/BOSC_Sternberg/B_analyses/B_BOSC/B_eBOSC_170703_full_wl6/';
% pn.tools = [pn.root, 'C_interrelations/T_tools/']; addpath(genpath(pn.tools));
% pn.XTask = [pn.root, 'B_extractIndices/B_data/'];
% pn.XRest = ['/Volumes/EEG/BOSC_SternRest/B_analyses/A_eBOSC/A_SternBRest_170703/B_extractIndices/B_data/'];
% pn.plotFolder = '/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F_Figures/';
% 
% Figure13.colorm = brewermap(4,'RdYlBu'); % use external function to create colormap
% 
% info.channels = 44:60;
% 
% load([pn.XRest, 'X15B_8to15_170703.mat']);
% 
% cond = {'EC1','EC2','EO1','EO2'};
% 
% for c = 1:numel(cond)
%     Figure13.Rest.abn_data(:,c) = nanmean(nanmean(X{1,1}.([cond{c},'e_abn'])(:,info.channels,:),3),2)';
%     % xAmp-BG
%     Figure13.Rest.oDiff(:,c) = nanmean(nanmean(X{1,1}.([cond{c},'o_amp'])(:,info.channels,:) - X{1,1}.([cond{c},'o_fitBG'])(:,info.channels,:),3),2)';
%     Figure13.Rest.eDiff(:,c) = nanmean(nanmean(X{1,1}.([cond{c},'e_amp'])(:,info.channels,:) - X{1,1}.([cond{c},'e_fitBG'])(:,info.channels,:),3),2)';
% end
% 
% %% get task data
% 
% pn.root = '/Volumes/EEG/BOSC_Sternberg/B_analyses/B_BOSC/B_eBOSC_170703_full_wl6/';
% pn.tools = [pn.root, 'C_interrelations/T_tools/']; addpath(genpath(pn.tools));
% pn.Xdata = [pn.root, 'B_extractIndices/B_data/'];
% 
% sessions = {'1'; '7'; '8'};
% indSession = 1;
% info.channels = 44:60;
% cond = {'L2','L4','L6'};
% 
% load([pn.XTask, 'X_8to15_170703_v3.mat']);
% 
% for c = 1:numel(cond)
%     Figure13.Task.abn_data(:,c) = nanmean(nanmean(X{1,1}.([cond{c},'e_abn'])(:,info.channels,:),3),2)';
%     % xAmp-BG
%     Figure13.Task.oDiff(:,c) = nanmean(nanmean(X{1,1}.([cond{c},'o_amp'])(:,info.channels,:) - X{1,1}.([cond{c},'o_fitBG'])(:,info.channels,:),3),2)';
%     Figure13.Task.aDiff(:,c) = nanmean(nanmean(X{1,1}.([cond{c},'a_amp'])(:,info.channels,:) - X{1,1}.([cond{c},'a_fitBG'])(:,info.channels,:),3),2)';
%     Figure13.Task.eDiff(:,c) = nanmean(nanmean(X{1,1}.([cond{c},'e_amp'])(:,info.channels,:) - X{1,1}.([cond{c},'e_fitBG'])(:,info.channels,:),3),2)';
% end
% 
% %% save Figure data
% 
% save('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/F13.mat', 'Figure13')

%% load Figure data

load('/Users/kosciessa/Desktop/eBOSC/figureData/F9C.mat', 'Figure13')

h = figure('units','normalized','position',[.1 .1 .25 .3]);
cla; hold on;
xData = 1-Figure13.Rest.abn_data(:,:);
yData = (Figure13.Rest.eDiff(:,:)-Figure13.Rest.oDiff(:,:))./Figure13.Rest.eDiff(:,:);
xData = squeeze(nanmean(xData(:,[1,3]),2)); yData = squeeze(nanmean(yData(:,[1,3]),2));
ls_1 = polyval(polyfit(xData(~isnan(xData) & ~isnan(yData)),yData(~isnan(xData) & ~isnan(yData)),1),0:1); ls_1 = plot(0:1, ls_1, 'Color', Figure13.colorm(1,:), 'LineWidth', 5);
scatter(xData(:),yData(:),150,Figure13.colorm(1,:), 'diamond', 'filled', 'CData', Figure13.colorm(1,:), 'MarkerEdgeColor', 'k', 'LineWidth',2);
[rl1,pl1] = corrcoef(xData,yData); pval1 = []; pval1 = convertPtoExponential(pl1(2));

xData = 1-Figure13.Rest.abn_data(:,:);
yData = (Figure13.Rest.eDiff(:,:)-Figure13.Rest.oDiff(:,:))./Figure13.Rest.eDiff(:,:);
xData = squeeze(nanmean(xData(:,[2,4]),2)); yData = squeeze(nanmean(yData(:,[2,4]),2));
ls_2 = polyval(polyfit(xData(~isnan(xData) & ~isnan(yData)),yData(~isnan(xData) & ~isnan(yData)),1),0:1); ls_2 = plot(0:1, ls_2, 'Color', Figure13.colorm(2,:), 'LineWidth', 5);
scatter(xData(:),yData(:),150,Figure13.colorm(2,:), 'diamond', 'filled', 'CData', Figure13.colorm(2,:), 'MarkerEdgeColor', 'k', 'LineWidth',2);
[rl2,pl2] = corrcoef(xData,yData); pval2 = []; pval2 = convertPtoExponential(pl2(2));

xData = 1-Figure13.Task.abn_data(:,:);
yData = (Figure13.Task.eDiff(:,:)-Figure13.Task.oDiff(:,:))./Figure13.Task.eDiff(:,:);
xData = squeeze(nanmean(xData(:,:),2)); yData = squeeze(nanmean(yData(:,:),2));
ls_3 = polyval(polyfit(xData(~isnan(xData) & ~isnan(yData)),yData(~isnan(xData) & ~isnan(yData)),1),0:1); ls_3 = plot(0:1, ls_3, 'Color', Figure13.colorm(4,:), 'LineWidth', 5);
scatter(xData(:),yData(:),150,Figure13.colorm(3,:), 'diamond', 'filled', 'CData', Figure13.colorm(4,:), 'MarkerEdgeColor', 'k', 'LineWidth',2);
[rl3,pl3] = corrcoef(xData,yData); pval3 = []; pval3 = convertPtoExponential(pl3(2));
hold on; plot([1 0],[1 0],'k--', 'LineWidth', 2)

% rearrange dots
childs = get(gca,'Children');
set(gca,'Children',[childs(2) childs(4) childs(6) childs(1) childs(3) childs(5) childs(7)])

lgd1 = legend([ls_1, ls_2, ls_3], {['S1 EC: r=', num2str(round(rl1(2),2)), ', p=',pval1{1}];...
    ['S1 EO: r=', num2str(round(rl2(2),2)), ', p=',pval2{1}];...
    ['S1 task: r=', num2str(round(rl3(2),2)), ', p=',pval3{1}]}, 'Location', 'NorthWest'); legend('boxoff');
xlabel('Arrhythmic duration (1-abundance)'); ylabel('Rhythm-specific amplitude gain')
title({'Arrhymic content linearly decreases';'traditional amplitude estimates'})
set(findall(gcf,'-property','FontSize'),'FontSize',20)
xlim([0 1]); ylim([0 1])

saveas(h, [pn.plotFolder, 'F13_v2'], 'fig');
saveas(h, [pn.plotFolder, 'F13_v2'], 'epsc');
saveas(h, [pn.plotFolder, 'F13_v2'], 'png');