%% Plot topographies of overall power, background and abundance + Spectra

% pn.root = '/Users/kosciessa/Desktop/mntTardis/Stern_WIP/B_eBOSC_180917_retention/';
% pn.tools = [pn.root, 'C_interrelations/T_tools/']; addpath(genpath(pn.tools));
% pn.Xdata = [pn.root, 'B_extractIndices/B_data/'];
% pn.plotFolder = [pn.root, 'X_AmpAbnTopographies/C_figures/'];
% pn.FieldTrip = [pn.root, 'A_eBOSC_180917_retention/T_tools/fieldtrip-20180227/'];
% addpath(pn.FieldTrip); ft_defaults;
% 
% %% load task data
% 
% load([pn.Xdata, 'X.mat']);
% 
% %% setup
% 
% sessions = {'1'; '7'; '8'};
% info.channels = 1:60;
% 
% cond = {'L2','L4','L6'};
% 
% Figure11.frequencies = {'Delta'; 'Theta'; 'Alpha'; 'Beta'; 'Gamma1'; 'Gamma2'};
% 
% for indSession = 1:3
% for indFreq = 1:numel(Figure11.frequencies)
% for c = 1:numel(cond)
%     % trial average
%     Figure11.o_data(indSession,indFreq,c,:,:) = nanmean(X{indSession,1}.(Figure11.frequencies{indFreq}).([cond{c},'o_amp'])(:,info.channels,:),3)';
%     Figure11.o_BGfit(indSession,indFreq,c,:,:) = nanmean(X{indSession,1}.(Figure11.frequencies{indFreq}).([cond{c},'o_fitBG'])(:,info.channels,:),3)';
%     Figure11.abn_data(indSession,indFreq,c,:,:) = nanmean(X{indSession,1}.(Figure11.frequencies{indFreq}).([cond{c},'e_abn'])(:,info.channels,:),3)';
% end % condition
% end % frequency
% end % session
% 
% dimord = 'sess_freq_cond_chan_sub';
% 
% 
% Figure11.cfg = [];
% Figure11.cfg.layout = 'elec1010.lay';
% Figure11.cfg.parameter = 'powspctrm';
% Figure11.cfg.colorbar = 'SouthOutside';
% Figure11.cfg.comment = 'no';
% Figure11.cfg.zlim = 'zeromax';
% 
% load('/Volumes/EEG/BOSC_SternRest/A_data/elec.mat');
% 
% Figure11.plotData = [];
% Figure11.plotData.label = elec.label(1:60,:); % {1 x N}
% Figure11.plotData.dimord = 'chan';
% 
% %% get spectra
% 
% for indSession = 1:3
% 
%     pn.bosc = [pn.root, 'A_eBOSC_180917_retention/B_data/B_BOSCout/eBOSC_S', sessions{indSession}, '/'];
% 
%     % N = 32 YA
%     IDs = {'3130';'3131';'3132';'3134';'3135';'3136';'3138';'3146';'3147';...
%         '3149';'3154';'3156';'3157';'3158';'3159';'3161';'3162';'3233';...
%         '3237';'3239';'3240';'3241';'3242';'3243';'3245';'3248';'3250';...
%         '3251';'3252';'3253';'3255';'3260'};
% 
%     fileNames = strcat(IDs, repmat('_bosc.mat',numel(IDs),1));
% 
%     IDs = cellfun(@str2num,IDs);
% 
%     cond = {'L2','L4','L6'};
%     condlab = {'load2', 'load4', 'load6'};
% 
%     for indID = 1:length(IDs)
%         % define ID
%         ID = num2str([IDs(indID,1)]); disp(ID);
% 
%         % load BGinfo
%         load([pn.bosc ID, '_bosc.mat'], 'BGinfo');
% 
%         Figure11.Thresholds.load2_pt(indSession,indID,:,:)        = BGinfo.all.load2_pt;
%         Figure11.Thresholds.load2_pv(indSession,indID,:,:)        = BGinfo.all.load2_pv;
%         Figure11.Thresholds.load2_bg_pow(indSession,indID,:,:)    = BGinfo.all.load2_bg_pow;
% 
%         Figure11.Thresholds.load4_pt(indSession,indID,:,:)        = BGinfo.all.load4_pt;
%         Figure11.Thresholds.load4_pv(indSession,indID,:,:)        = BGinfo.all.load4_pv;
%         Figure11.Thresholds.load4_bg_pow(indSession,indID,:,:)    = BGinfo.all.load4_bg_pow;
% 
%         Figure11.Thresholds.load6_pt(indSession,indID,:,:)        = BGinfo.all.load6_pt;
%         Figure11.Thresholds.load6_pv(indSession,indID,:,:)        = BGinfo.all.load6_pv;
%         Figure11.Thresholds.load6_bg_pow(indSession,indID,:,:)    = BGinfo.all.load6_bg_pow;
%     end
%     
% end
% 
% % average across sessions and loads
% Figure11.Thresholds.allLoads_bg_pow = cat(1, Figure11.Thresholds.load2_bg_pow, Figure11.Thresholds.load4_bg_pow, Figure11.Thresholds.load6_bg_pow);
% Figure11.Thresholds.allLoads_bg_pow = squeeze(nanmean(Figure11.Thresholds.allLoads_bg_pow,1));
% 
% Figure11.freq = 2.^[0:.125:6];
% 
% %% save Figure data
% 
% save('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/F11.mat', 'Figure11')

%% load Figure data

addpath('/Volumes/EEG/BOSC_Sternberg/T_tools/fieldtrip-20180227'); ft_defaults;

load('/Users/kosciessa/Desktop/eBOSC/figureData/F8.mat', 'Figure11')

%% plot average across loads & sessions: power threshold

% average: cond, sess, sub for each Figure11.freq

h = figure('units','normalized','position',[.1 .1 .85 .4]);

indFreq = 3; % 8-15 Hz Alpha

subplot(1,4,1);
    Figure11.cfg.highlight = 'off';
    Figure11.plotData.powspctrm = squeeze(nanmean(nanmean(nanmean(Figure11.o_data(1:3,indFreq,:,:,:),3),1),5));
    ft_topoplotER(Figure11.cfg,Figure11.plotData); 
    title(['Overall ',Figure11.frequencies{indFreq},' amplitude']);
subplot(1,4,2);
    Figure11.cfg.highlight = 'on';
    Figure11.cfg.highlightcolor = [1 1 1];
    Figure11.cfg.highlightchannel = 'FCz';
    Figure11.cfg.highlightsize = 10;
    Figure11.plotData.powspctrm = squeeze(nanmean(nanmean(nanmean(Figure11.o_BGfit(1:3,indFreq,:,:,:),3),1),5));
    ft_topoplotER(Figure11.cfg,Figure11.plotData)
    title([Figure11.frequencies{indFreq},' background amplitude']);
subplot(1,4,3);
    Figure11.cfg.highlight = 'on';
    Figure11.cfg.highlightsize = 10;
    Figure11.cfg.highlightcolor = [1 1 1];
    Figure11.cfg.highlightchannel = 44:60;
    Figure11.plotData.powspctrm = squeeze(nanmean(nanmean(nanmean(Figure11.abn_data(1:3,indFreq,:,:,:),3),1),5));    
    ft_topoplotER(Figure11.cfg,Figure11.plotData)
    title([Figure11.frequencies{indFreq},' abundance']);
 
% add colorbrewer
addpath('/Volumes/LNDG/Projects/StateSwitch/dynamic/data/eeg/rest/B_analyses/A_MSE_CSD_multiVariant/T_tools/brewermap')
cBrew = brewermap(500,'RdBu');
cBrew = flipud(cBrew);
colormap(cBrew)

%% plot spectra separately for frontal/posterior channels

pn.shadedError = ['/Volumes/LNDG/Projects/StateSwitch/dynamic/data/eeg/rest/B_analyses/A_MSE_CSD_multiVariant/T_tools/shadedErrorBar']; addpath(pn.shadedError);

subplot(1,4,4); cla; hold on;
curAverage = nanmean(Figure11.Thresholds.allLoads_bg_pow(:,21,:),2);
    standError = nanstd(curAverage,1)./sqrt(size(curAverage,1));
    l1 = shadedErrorBar([],nanmean(curAverage,1),standError, 'lineprops', {'Color',[.3 .5 .8],'linewidth', 4}, 'patchSaturation', .05);
curAverage = nanmean(Figure11.Thresholds.allLoads_bg_pow(:,44:60,:),2);
    standError = nanstd(curAverage,1)./sqrt(size(curAverage,1));
    l2 = shadedErrorBar([],nanmean(curAverage,1),standError, 'lineprops', {'r','linewidth', 4}, 'patchSaturation', .05);

% plot(squeeze(nanmean(nanmean(Figure11.Thresholds.allLoads_bg_pow(:,21,:),2),1)), 'LineWidth', 4);
% plot(squeeze(nanmean(nanmean(Figure11.Thresholds.allLoads_bg_pow(:,44:60,:),2),1)), 'LineWidth', 4);
xlim([1, numel(Figure11.freq)])
xticks([1:10:numel(Figure11.freq)]); xticklabels(round(Figure11.freq(xticks),1)); xlabel('Frequency [Hz]')
legend([l1.mainLine, l2.mainLine], {'FCz average spectrum'; 'Posterior average spectrum'}); legend('boxoff')
title({'Spectra for labeled channels'; ''})
set(findall(gcf,'-property','FontSize'),'FontSize',20)

h_sup = suptitle('eBOSC differentiates between frontal 1/f slope and posterior rhythmicity during WM retention');
set(h_sup, 'FontWeight','bold', 'FontSize', 24);

pn.plotFolder = '/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/Figures_resubmit1/F1_Figures/';
figureName = 'F8';

saveas(h, [pn.plotFolder, figureName], 'fig');
saveas(h, [pn.plotFolder, figureName], 'epsc');
saveas(h, [pn.plotFolder, figureName], 'png');

