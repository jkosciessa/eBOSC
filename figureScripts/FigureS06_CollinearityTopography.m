%% Plot collinearity topography between overall power and abundance

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
% cond = {'L2','L4','L6'};
% FigureS5.frequencies = {'Delta'; 'Theta'; 'Alpha'; 'Beta'; 'Gamma1'; 'Gamma2'};
% 
% for indSession = 1:3
% for indFreq = 1:numel(FigureS5.frequencies)
% for c = 1:numel(cond)
%     o_data(indSession,indFreq,c,:,:) = nanmean(X{indSession,1}.(FigureS5.frequencies{indFreq}).([cond{c},'o_amp'])(:,1:60,:),3)';
%     o_data_ST(indSession,indFreq,c,:,:,:) = X{indSession,1}.(FigureS5.frequencies{indFreq}).([cond{c},'o_amp'])(:,1:60,:);
%     abn_data(indSession,indFreq,c,:,:) = nanmean(X{indSession,1}.(FigureS5.frequencies{indFreq}).([cond{c},'e_abn'])(:,1:60,:),3)';
%     abn_data_ST(indSession,indFreq,c,:,:,:) = X{indSession,1}.(FigureS5.frequencies{indFreq}).([cond{c},'e_abn'])(:,1:60,:);
% end % condition
% end % frequency
% end % session
% FigureS5.dimord = 'sess_freq_cond_chan_sub';
% 
% FigureS5.cfg = [];
% FigureS5.cfg.layout = 'elec1010.lay';
% FigureS5.cfg.comment = 'no';
% 
% load('/Volumes/EEG/BOSC_SternRest/A_data/elec.mat');
% 
% FigureS5.plotData = [];
% FigureS5.plotData.label = elec.label(1:60,:); % {1 x N}
% FigureS5.plotData.dimord = 'chan';
% 
% % between-subject collinearity
% 
% for indFreq = 1:4
%     o_data_BySub = squeeze(nanmean(nanmean(o_data(:,indFreq,:,:,:),3),1));
%     abn_data_BySub = squeeze(nanmean(nanmean(abn_data(:,indFreq,:,:,:),3),1));
%     for indChan = 1:60
%         [r] = corrcoef(o_data_BySub(indChan,:), abn_data_BySub(indChan,:));
%         FigureS5.BetweenCollinearityByFreq{indFreq}(indChan,1) = r(2);
%     end
% end
% 
% % within-subject collinearity
% 
% for indFreq = 1:4
%     for indSub = 1:32
%         for indChan = 1:60
%             tmp_o_data_BySub = o_data_ST(:,indFreq,:,indSub,indChan,:);
%             tmp_o_data_BySub = tmp_o_data_BySub(~isnan(tmp_o_data_BySub));
%             tmp_abn_data_BySub = abn_data_ST(:,indFreq,:,indSub,indChan,:);
%             tmp_abn_data_BySub = tmp_abn_data_BySub(~isnan(tmp_abn_data_BySub));
%             [r] = corrcoef(tmp_o_data_BySub, tmp_abn_data_BySub);
%             FigureS5.WithinCollinearityByFreq{indFreq}(indChan,indSub) = r(2);
%         end
%     end
%     FigureS5.WithinCollinearityByFreq{indFreq} = squeeze(nanmean(FigureS5.WithinCollinearityByFreq{indFreq},2));
% end
% 
% %% save Figure data
% 
% save('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/S5.mat', 'FigureS5')

%% load Figure data

load('/Users/kosciessa/Desktop/eBOSC/figureData/S6.mat', 'FigureS5')

%% plot topography of between-subject collinearity

h = figure('units','normalized','position',[.1 .1 .7 .5]);

FigureS5.cfg.zlim = [-1 1];
FigureS5.cfg.parameter = 'BetweenCollinearity';

for indFreq = 1:4
    subplot(2,4,indFreq);
    FigureS5.plotData.BetweenCollinearity = FigureS5.BetweenCollinearityByFreq{indFreq};
    ft_topoplotER(FigureS5.cfg,FigureS5.plotData); colorbar('location', 'SouthOutside');
    title({'Between-subject:'; FigureS5.frequencies{indFreq}});
end

%% plot topography of within-subject collinearity

FigureS5.cfg.zlim = [-.6 .6];
FigureS5.cfg.parameter = 'WithinCollinearity';

for indFreq = 1:4
    subplot(2,4,4+indFreq);
    FigureS5.plotData.WithinCollinearity = FigureS5.WithinCollinearityByFreq{indFreq};
    ft_topoplotER(FigureS5.cfg,FigureS5.plotData); colorbar('location', 'SouthOutside');
    title({'Within-subject:'; FigureS5.frequencies{indFreq}});
end

suptitle('Collinearity between overall amplitude & abundance')

pn.plotFolder = '/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F_Figures/';
figureName = 'S5';

saveas(h, [pn.plotFolder, figureName], 'fig');
saveas(h, [pn.plotFolder, figureName], 'epsc');
saveas(h, [pn.plotFolder, figureName], 'png');