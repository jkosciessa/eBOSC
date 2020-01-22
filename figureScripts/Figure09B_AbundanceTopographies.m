%% Plot abundance topographies for multiple frequency intervals

% pn.root = '/Volumes/EEG/BOSC_Sternberg/B_analyses/B_BOSC/B_eBOSC_180917_retention/';
% pn.abnData = [pn.root, 'X_AmpAbnTopographies/B_data/'];
% pn.FieldTrip = [pn.root, 'A_eBOSC_180917_retention/T_tools/fieldtrip-20180227/'];
% addpath(pn.FieldTrip); ft_defaults;
% 
% pn.plotFolder = [pn.root, 'X_AmpAbnTopographies/C_figures/'];
% 
% %% load abundance-by-channel data
% 
% load([pn.abnData, 'Abundance@allChannnels.mat'], 'Abundance');
% 
% %% load channel layout
% 
% load('/Volumes/EEG/BOSC_Sternberg/A_data/JK_periRet_Ext_corr_1/3260-01-ArtCorr.mat'); % needs the original data structure
% layoutInfo.elec = data.elec; clear data;
% layoutInfo.elec.chanpos = layoutInfo.elec.chanpos(1:60,:);
% layoutInfo.elec.elecpos = layoutInfo.elec.elecpos(1:60,:);
% layoutInfo.elec.label = layoutInfo.elec.label(1:60,:);
% 
% %% collapse across sessions
% 
% conditions = {'load2', 'load4', 'load6'};
% Figure12B.FrequencyCategories = {'Delta'; 'Theta'; 'Alpha'; 'Beta'; 'Gamma1'; 'Gamma2'};
% Figure12B.FrequencyLabels = {{'Delta';'(2.5 - 4 Hz)'}; {'Theta';'(4 - 8 Hz)'}; {'Alpha';'(8 - 15 Hz)'};...
%     {'Beta';'(15 - 25 Hz)'}; {'Gamma 1';'(25 - 40 Hz)'}; {'Gamma 2';'(40 - 64 Hz)'}};
% 
% for indCond = 1:numel(conditions)
%     for indFreq = 1:numel(Figure12B.FrequencyCategories)
%         curAbn = Abundance.(Figure12B.FrequencyCategories{indFreq});
%         AbundanceMerged.(Figure12B.FrequencyCategories{indFreq}).(conditions{indCond}) = ...
%             cat(4, curAbn{1,1}.(conditions{indCond}), curAbn{2,1}.(conditions{indCond}), ...
%             curAbn{3,1}.(conditions{indCond}));
%         AbundanceMerged.(Figure12B.FrequencyCategories{indFreq}).(conditions{indCond}) = ...
%             nanmean(AbundanceMerged.(Figure12B.FrequencyCategories{indFreq}).(conditions{indCond}),4);
%         % average across trials
%         AbundanceMerged.(Figure12B.FrequencyCategories{indFreq}).(conditions{indCond}) = ...
%             squeeze(nanmean(AbundanceMerged.(Figure12B.FrequencyCategories{indFreq}).(conditions{indCond}),2));
%     end
% end
% 
% %% plot aggregated across conditions
% 
% Figure12B.cfg = [];
% Figure12B.cfg.layout = 'elec1010.lay';
% Figure12B.cfg.parameter = 'powspctrm';
% Figure12B.cfg.colorbar = 'South';
% Figure12B.cfg.comment = 'no';
% Figure12B.cfg.zlim = 'zeromax';
% 
% Figure12B.plotData = [];
% Figure12B.plotData.label = layoutInfo.elec.label'; % {1 x N}
% Figure12B.plotData.dimord = 'chan';
% 
% for indFreq = 1:numel(Figure12B.FrequencyCategories) % average across loads
%     Figure12B.AbundanceByFreq{indFreq} = squeeze(nanmean(cat(3, AbundanceMerged.(Figure12B.FrequencyCategories{indFreq}).load2, ...
%         AbundanceMerged.(Figure12B.FrequencyCategories{indFreq}).load4,...
%         AbundanceMerged.(Figure12B.FrequencyCategories{indFreq}).load6),3));
% end
% 
% %% save Figure data
% 
% save('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/F12B.mat', 'Figure12B')

%% load Figure data

load('/Users/kosciessa/Desktop/eBOSC/figureData/F9B.mat', 'Figure12B')

%% plot Abundance topographies

h = figure('units','normalized','position',[.1 .1 .7 .3]);
for indFreq = 1:numel(Figure12B.FrequencyCategories)
    subplot(1,6,indFreq)
    Figure12B.cfg.colorbar = 'SouthOutside'
    Figure12B.plotData.powspctrm = squeeze(nanmean(Figure12B.AbundanceByFreq{indFreq},1))';
    ft_topoplotER(Figure12B.cfg,Figure12B.plotData)
    title(Figure12B.FrequencyLabels{indFreq})
    cb = colorbar('location', Figure12B.cfg.colorbar);
    set(get(cb,'XLabel'),'String','abundance');
end
set(findall(gcf,'-property','FontSize'),'FontSize',25)

% change colormap
addpath('/Volumes/EEG/BOSC_Sternberg/T_tools/brewermap')
cBrew = brewermap(500,'RdBu');
cBrew = flipud(cBrew);
colormap(cBrew)

pn.plotFolder = '/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F_Figures/';
figureName = 'F12B';

saveas(h, [pn.plotFolder, figureName], 'fig');
saveas(h, [pn.plotFolder, figureName], 'epsc');
saveas(h, [pn.plotFolder, figureName], 'png');
