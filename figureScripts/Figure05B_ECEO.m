% - used to be /Volumes/EEG/BOSC_SternRest/B_scripts_JQK/B_processing_170704/E_JQKplotAmpAbn_170816B.mlx

% pn.root     = '/Volumes/EEG/BOSC_SternRest/B_analyses/A_eBOSC/A_SternBRest_170703/';
% pn.tools    = [pn.root, 'C_interrelations/T_tools/']; addpath(genpath(pn.tools));
% pn.data     = [pn.root, 'B_extractIndices/B_data/'];
% 
% %% load data
% 
% load([pn.data, 'X15B_8to15_170703.mat']);
% 
% %% extract relevant metrics
% 
% info.channels = 44:60;
% cond = {'EC1','EC2','EO1','EO2'};
% for c = 1:numel(cond)
%     abn_data(:,c) = nanmean(nanmean(X{1,1}.([cond{c},'e_abn'])(:,info.channels,:),3),2)';
%     % xAmp-BG
%     eDiff(:,c) = nanmean(nanmean(X{1,1}.([cond{c},'e_amp'])(:,info.channels,:) - X{1,1}.([cond{c},'e_fitBG'])(:,info.channels,:),3),2)';
% end;
% 
% Figure5B.colorm = brewermap(4,'RdYlBu'); % use external function to create colormap
% 
% Figure5B.EO1e_abn = abn_data(:,3); Figure5B.EO1e_amp = eDiff(:,3);
% Figure5B.EC1e_abn = abn_data(:,1); Figure5B.EC1e_amp = eDiff(:,1);
% 
% Figure5B.EO2e_abn = abn_data(:,4); Figure5B.EO2e_amp = eDiff(:,4);
% Figure5B.EC2e_abn = abn_data(:,2); Figure5B.EC2e_amp = eDiff(:,2);
% 
% %% save Figure data
% 
% save('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/F5B.mat', 'Figure5B')

%% load Figure data

load('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/F5B.mat', 'Figure5B')

%% Figure: Eye closure state change

h = figure('units','normalized','position',[.1 .1 .3 .7]);
title('Parameter change upon eye closure (individual subjects)')
for j = 1:32 % EO1-EC1
    hold on; h1 = arrow([Figure5B.EO1e_abn(j) Figure5B.EO1e_amp(j)],[Figure5B.EC1e_abn(j) Figure5B.EC1e_amp(j)], 5, 'Color',Figure5B.colorm(1,:));
end; clear j
hold on; arrow([mean(Figure5B.EO1e_abn) mean(Figure5B.EO1e_amp)],[mean(Figure5B.EC1e_abn) mean(Figure5B.EC1e_amp)], 5, 'Color',[.5 0 0], 'LineWidth', 3);
for j = 1:32 % EO2-EC2
    hold on; h2 = arrow([Figure5B.EO2e_abn(j) Figure5B.EO2e_amp(j)],[Figure5B.EC2e_abn(j) Figure5B.EC2e_amp(j)], 5, 'Color',Figure5B.colorm(4,:));
end; clear j
hold on; arrow([mean(Figure5B.EO2e_abn) mean(Figure5B.EO2e_amp)],[mean(Figure5B.EC2e_abn) mean(Figure5B.EC2e_amp)], 5, 'Color',[0 0 .5], 'LineWidth', 3);
xlabel('Abundance'); ylabel('Rhythmic Amplitude excl. BG [?V]');
legend([h1, h2], {'continuous', 'interleaved'}, 'location', 'northwest');
legend('boxoff')
set(findall(gcf,'-property','FontSize'),'FontSize',15)

pn.plotFolder = '/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F_Figures/';

saveas(h, [pn.plotFolder, 'F5B'], 'fig');
saveas(h, [pn.plotFolder, 'F5B'], 'epsc');
