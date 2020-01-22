%% Plot overall amplitude, background amplitude and abundance topographies for multiple frequencies

% pn.root = '/Users/kosciessa/Desktop/mntTardis/Stern_WIP/B_eBOSC_180917_retention/';
% pn.tools = [pn.root, 'C_interrelations/T_tools/']; addpath(genpath(pn.tools));
% pn.Xdata = [pn.root, 'B_extractIndices/B_data/'];
% 
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
% FigureS4.frequencies = {'Delta'; 'Theta'; 'Alpha'; 'Beta'; 'Gamma1'; 'Gamma2'};
% 
% for indSession = 1:3
% for indFreq = 1:numel(FigureS4.frequencies)
% for c = 1:numel(cond)
%     % trial average
%     o_data(indSession,indFreq,c,:,:) = nanmean(X{indSession,1}.(FigureS4.frequencies{indFreq}).([cond{c},'o_amp'])(:,info.channels,:),3)';
%     o_BGfit(indSession,indFreq,c,:,:) = nanmean(X{indSession,1}.(FigureS4.frequencies{indFreq}).([cond{c},'o_fitBG'])(:,info.channels,:),3)';
%     abn_data(indSession,indFreq,c,:,:) = nanmean(X{indSession,1}.(FigureS4.frequencies{indFreq}).([cond{c},'e_abn'])(:,info.channels,:),3)';
% end % condition
% end % frequency
% end % session
% 
% dimord = 'sess_freq_cond_chan_sub';
% 
% FigureS4.cfg = [];
% FigureS4.cfg.layout = 'elec1010.lay';
% FigureS4.cfg.parameter = 'powspctrm';
% FigureS4.cfg.colorbar = 'East';
% FigureS4.cfg.comment = 'no';
% FigureS4.cfg.zlim = 'zeromax';
% 
% load('/Volumes/EEG/BOSC_SternRest/A_data/elec.mat');
% 
% FigureS4.plotData = [];
% FigureS4.plotData.label = elec.label(1:60,:); % {1 x N}
% FigureS4.plotData.dimord = 'chan';
% 
% % average across sessions, loads & subjects
% for indFreq = 1:numel(FigureS4.frequencies)
%     FigureS4.oDataMergedbyFreq{indFreq} = ...
%         squeeze(nanmean(nanmean(nanmean(o_data(:,indFreq,:,:,:),3),1),5));
%     FigureS4.oBGfitMergedbyFreq{indFreq} = ...
%         squeeze(nanmean(nanmean(nanmean(o_BGfit(:,indFreq,:,:,:),3),1),5));
%     FigureS4.abnDataMergedbyFreq{indFreq} = ...
%         squeeze(nanmean(nanmean(nanmean(abn_data(:,indFreq,:,:,:),3),1),5));
% end
% 
% %% save Figure data
% 
% save('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/S4.mat', 'FigureS4')

%% load Figure data

load('/Users/kosciessa/Desktop/eBOSC/figureData/S8.mat', 'FigureS4')

%% plot overall amplitude, background amplitude and abundance

h = figure('units','normalized','position',[.1 .1 .7 .9]);
for indFreq = 1:numel(FigureS4.frequencies)
    subplot(6,3,(indFreq-1)*3+1);
        FigureS4.plotData.powspctrm = FigureS4.oDataMergedbyFreq{indFreq};
        ft_topoplotER(FigureS4.cfg,FigureS4.plotData); 
        title(['Overall ',FigureS4.frequencies{indFreq},' amplitude']);
    subplot(6,3,(indFreq-1)*3+2);
        FigureS4.plotData.powspctrm = FigureS4.oBGfitMergedbyFreq{indFreq};
        ft_topoplotER(FigureS4.cfg,FigureS4.plotData)
        title([FigureS4.frequencies{indFreq},' background amplitude']);
    subplot(6,3,(indFreq-1)*3+3);
        FigureS4.plotData.powspctrm = FigureS4.abnDataMergedbyFreq{indFreq};
        ft_topoplotER(FigureS4.cfg,FigureS4.plotData)
        title([FigureS4.frequencies{indFreq},' abundance']);
end

pn.plotFolder = '/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F_Figures/';
figureName = 'S4';

saveas(h, [pn.plotFolder, figureName], 'fig');
saveas(h, [pn.plotFolder, figureName], 'epsc');
saveas(h, [pn.plotFolder, figureName], 'png');
