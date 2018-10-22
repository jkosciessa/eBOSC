%% Plot topographies for resting state + state modulation

% original: /Volumes/EEG/BOSC_SternRest/B_scripts_JQK/B_processing_170704/H_plotAbundanceTopography_170816B2.mlx

% clear all; clc; restoredefaultpath;
% 
% addpath('/Volumes/EEG/ConMemEEGTools/fieldtrip-20180227/')
% 
% load('/Volumes/EEG/BOSC_SternRest/A_data/elec.mat');
% load('/Volumes/EEG/BOSC_SternRest/B_analyses/A_eBOSC/A_SternBRest_170703/B_extractIndices/B_data/X15B_8to15_170703.mat');
% 
% %% create Figure5A.plotData structure
% 
% % average across trials and subjects
% 
% Figure5A.plotData = [];
% 
% Figure5A.plotData.AlphaAbn_EO = squeeze(nanmean(nanmean(cat(3, X{1,1}.EO1e_abn, X{1,1}.EO2e_abn),3),1))';
% Figure5A.plotData.AlphaAbn_EC = squeeze(nanmean(nanmean(cat(3, X{1,1}.EC1e_abn, X{1,1}.EC2e_abn),3),1))';
% 
% Figure5A.plotData.AlphaAmp_EO = squeeze(nanmean(nanmean(cat(3, X{1,1}.EO1e_amp_BGdiff, X{1,1}.EO2e_amp_BGdiff),3),1))';
% Figure5A.plotData.AlphaAmp_EC = squeeze(nanmean(nanmean(cat(3, X{1,1}.EC1e_amp_BGdiff, X{1,1}.EC2e_amp_BGdiff),3),1))';
% 
% Figure5A.plotData.AlphaAmp_a_EO = squeeze(nanmean(nanmean(cat(3,X{1,1}.EO1a_amp_BGdiff,X{1,1}.EO2a_amp_BGdiff),3),1))';
% Figure5A.plotData.AlphaAmp_a_EC = squeeze(nanmean(nanmean(cat(3,X{1,1}.EC1a_amp_BGdiff,X{1,1}.EC2a_amp_BGdiff),3),1))';
% 
% Figure5A.plotData.AlphaAmp_na_EO = squeeze(nanmean(nanmean(cat(3,X{1,1}.EO1na_amp_BGdiff,X{1,1}.EO2na_amp_BGdiff),3),1))';
% Figure5A.plotData.AlphaAmp_na_EC = squeeze(nanmean(nanmean(cat(3,X{1,1}.EC1na_amp_BGdiff,X{1,1}.EC2na_amp_BGdiff),3),1))';
% 
% Figure5A.plotData.label = elec.label(1:60,:); % {1 x N}
% Figure5A.plotData.dimord = 'chan';
% 
% %% create Figure5A.cfg
% 
% Figure5A.cfg = [];
% Figure5A.cfg.layout = 'elec1010.lay';
% Figure5A.cfg.colorbar = 'no';
% Figure5A.cfg.comment = 'no';
% 
% save('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/F5A.mat', 'Figure5A')

%% load Figure data

load('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/F5A.mat', 'Figure5A')

%% plot aggregated across conditions

h = figure('units','normalized','position',[.1 .1 .4 .7]);

conds = {'EO';'EC'};
condLabels = {'Eyes Open';'Eyes Closed'};

for indCondition = 1:numel(conds)
    Figure5A.cfg.parameter = ['AlphaAbn_',conds{indCondition}];
    subplot(3,numel(conds),indCondition);
    if indCondition == 1
        Figure5A.cfg.zlim = [0 .4];
    else
        Figure5A.cfg.zlim = [0 .8];
    end
    ft_topoplotER(Figure5A.cfg,Figure5A.plotData); colorbar('location', 'EastOutside');
    if indCondition == 1
        text(-0.4,0.5,{'Abundance'},'units','normalized', 'FontSize', 13, 'FontWeight', 'bold');
    end;
    text(0.3,1.15,condLabels{indCondition},'units','normalized', 'FontSize', 13, 'FontWeight', 'bold');
end;

for indCondition = 1:numel(conds)
    Figure5A.cfg.parameter = ['AlphaAmp_',conds{indCondition}];
    subplot(3,numel(conds),numel(conds)+indCondition);
    if indCondition == 1
        Figure5A.cfg.zlim = [0 500];
    else
        Figure5A.cfg.zlim = [0 1000];
    end
    ft_topoplotER(Figure5A.cfg,Figure5A.plotData); colorbar('location', 'EastOutside');
    if indCondition == 1
        text(-0.4,0.5,{'Rhythmic','amplitude','(excl. BG)'},'units','normalized', 'FontSize', 13, 'FontWeight', 'bold');
    end;
end;

for indCondition = 1:numel(conds)
    Figure5A.cfg.parameter = ['AlphaAmp_na_',conds{indCondition}];
    subplot(3,numel(conds),2*numel(conds)+indCondition);
    if indCondition == 1
        Figure5A.cfg.zlim = [0 500];
    else
        Figure5A.cfg.zlim = [0 1000];
    end
    ft_topoplotER(Figure5A.cfg,Figure5A.plotData); colorbar('location', 'EastOutside');
    if indCondition == 1
        text(-0.4,0.5,{'Arrhythmic';'amplitude'},'units','normalized', 'FontSize', 13, 'FontWeight', 'bold');
    end;
end;

set(findall(gcf,'-property','FontSize'),'FontSize',15)

pn.plotFolder = '/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F_Figures/';

saveas(h, [pn.plotFolder, 'F5A'], 'fig');
saveas(h, [pn.plotFolder, 'F5A'], 'epsc');
