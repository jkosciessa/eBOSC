% Calculate amplitude-abundance link across load levels

% based originally on /Volumes/EEG/BOSC/BOSC_Sternberg/scripts/B_BOSC/C_JQKplotAmpAbn_170703.mlx

% 180920 | use updated X structure with correct alpha range

% %% set paths
% 
% pn.root = '/Volumes/EEG/BOSC/BOSC_Sternberg/B_analyses/B_BOSC/B_eBOSC_170703_full_wl6/';
% pn.tools = [pn.root, 'C_interrelations/T_tools/']; addpath(genpath(pn.tools));
% pn.Xdata = [pn.root, 'B_extractIndices/B_data/'];
% 
% %% setup
% 
% sessions = {'1'; '7'; '8'};
% info.channels = 44:60;
% Figure8.colormBetweenSub = brewermap(4,'RdYlBu'); % use external function to create colormap
% Figure8.colormWithinSub = brewermap(32,'RdYlBu'); % use external function to create colormap
% 
% %% collect data
% 
% load([pn.Xdata, 'X_8to15_170703_v3.mat']);
% 
% % IMPORTANT: Note that variables of interest have to be computed on the
% % trial-level first. Working from averages will give different results as
% % there are missings, which will alter the amount of variables.
% 
% % Current approach: average within each subject, channel & trial produces
% % the matrices that averages will be computed on below.
% 
% Figure8.cond = {'L2','L4','L6'};
% 
% for indSession = 1
%     for c = 1:numel(Figure8.cond)
%         % single-trial amplitudes (avg. across channels); note there are no single trials for o & a
%         Figure8.e_amp_singleTrial(indSession,:,c,:) = nanmean(X{indSession,1}.([Figure8.cond{c},'e_amp'])(:,info.channels,:),2);
%         Figure8.e_BG_singleTrial(indSession,:,c,:) = nanmean(X{indSession,1}.([Figure8.cond{c},'e_fitBG'])(:,info.channels,:),2);
%         % single-trial abundance (avg. across channels)
%         Figure8.e_abn_singleTrial(indSession,:,c,:) = nanmean(X{indSession,1}.([Figure8.cond{c},'e_abn'])(:,info.channels,:),2);
%         Figure8.eDiffRel_singleTrial(indSession,:,c,:) = nanmean((X{indSession,1}.([Figure8.cond{c},'e_amp'])(:,info.channels,:)-X{indSession,1}.([Figure8.cond{c},'e_fitBG'])(:,info.channels,:))./X{indSession,1}.([Figure8.cond{c},'e_fitBG'])(:,info.channels,:),2);
%     end
% end
% 
% %% save Figure data
% 
% save('/Volumes/EEG/BOSC/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/F8.mat', 'Figure8')

%% load Figure data

addpath('/Volumes/EEG/BOSC/BOSC_SternRest/X_documentation/B_2018_Manuscript/F4_FigureTools')

load('/Users/kosciessa/Desktop/eBOSC/figureData/F7.mat', 'Figure8')

%% Plot v2

h = figure('DefaultAxesFontSize',12, 'units','normalized','position',[.1 .1 .7 .3]);
subplot(1,3,1); cla;hold on;
    cmap = lines(32);
    cmap2 = brighten(cmap,.7);
    for indSub = 1:32
        xData = Figure8.e_abn_singleTrial(1,indSub,:,:);
        yData = Figure8.eDiffRel_singleTrial(1,indSub,:,:);
        ls{indSub} = polyval(polyfit(xData(~isnan(xData) & ~isnan(yData)),yData(~isnan(xData) & ~isnan(yData)),1),0:1); plot(0:1, ls{indSub}, 'Color', [.8 .8 .8], 'LineWidth', 1);
        scatter(xData(:),yData(:),75,Figure8.colormWithinSub(indSub,:),'filled', 'MarkerFaceColor', cmap2(indSub,:));
    end
    % add between-subject mean
    xData = squeeze(nanmean(nanmean(Figure8.e_abn_singleTrial(1,:,:,:),3),4));
    yData = squeeze(nanmean(nanmean(Figure8.eDiffRel_singleTrial(1,:,:,:),3),4));
    ls_all = polyval(polyfit(xData(~isnan(xData) & ~isnan(yData)),yData(~isnan(xData) & ~isnan(yData)),1),0:1); ls_all = plot(0:1, ls_all, 'Color', [0 0 0], 'LineWidth', 3);
    scatter(xData(:),yData(:),250,Figure8.colormWithinSub(indSub,:), 'diamond', 'filled', 'CData', cmap(:,:), 'MarkerEdgeColor', 'k', 'LineWidth',2);
    [rl1,pl1] = corrcoef(xData,yData); pval1 = []; pval1 = convertPtoExponential(pl1(2));
    lgd1 = legend([ls_all], {['r=', num2str(round(rl1(2),2)), ', p=',pval1{1}]}, 'Location', 'NorthWest'); legend('boxoff');
    title({'Rhythmic SNR & abundance relate ';'between and within subjects'}, 'FontSize', 20);
    xlabel('Rhythmic abundance'); ylabel({'Rhythmic SNR'});
    xlim([0 1]); ylim([0 12])
subplot(1,3,2); cla;hold on;
    cmap = lines(32);
    cmap2 = brighten(cmap,.7);
    for indSub = 1:32
        xData = Figure8.e_abn_singleTrial(1,indSub,:,:);
        yData = Figure8.e_BG_singleTrial(1,indSub,:,:);
        ls{indSub} = polyval(polyfit(xData(~isnan(xData) & ~isnan(yData)),yData(~isnan(xData) & ~isnan(yData)),1),0:1); plot(0:1, ls{indSub}, 'Color', [.8 .8 .8], 'LineWidth', 1);
        scatter(xData(:),yData(:),75,Figure8.colormWithinSub(indSub,:),'filled', 'MarkerFaceColor', cmap2(indSub,:));
    end
    % add between-subject mean
    xData = squeeze(nanmean(nanmean(Figure8.e_abn_singleTrial(1,[1:3,5:32],:,:),3),4));
    yData = squeeze(nanmean(nanmean(Figure8.e_BG_singleTrial(1,[1:3,5:32],:,:),3),4));
    xDataOutlier = squeeze(nanmean(nanmean(Figure8.e_abn_singleTrial(1,[4],:,:),3),4));
    yDataOutlier = squeeze(nanmean(nanmean(Figure8.e_BG_singleTrial(1,[4],:,:),3),4));
    ls_all = polyval(polyfit(xData(~isnan(xData) & ~isnan(yData)),yData(~isnan(xData) & ~isnan(yData)),1),0:1); ls_all = plot(0:1, ls_all, 'Color', [0 0 0], 'LineWidth', 3);
    scatter(xData(:),yData(:),250,Figure8.colormWithinSub(indSub,:), 'diamond', 'filled', 'CData', cmap([1:3,5:32],:), 'MarkerEdgeColor', 'k', 'LineWidth',2);
    scatter(xDataOutlier,yDataOutlier,250,Figure8.colormWithinSub(indSub,:), 'diamond', 'filled', 'CData', cmap(4,:), 'MarkerEdgeColor', [.9 .9 .9],'LineWidth',2);
    [rl1,pl1] = corrcoef(xData,yData); pval1 = []; pval1 = convertPtoExponential(pl1(2));
    lgd1 = legend([ls_all], {['r=', num2str(round(rl1(2),2)), ', p=',pval1{1}]}, 'Location', 'NorthWest'); legend('boxoff');
    title({'Background amplitude does not';'consistently relate to abundance'}, 'FontSize', 20);
    xlabel('Rhythmic abundance'); ylabel({'Background amplitude'});
    xlim([0 1]); ylim([75 220])
set(findall(gcf,'-property','FontSize'),'FontSize',22)

% add histogram for individual amplitude-abundance correlations

handaxes1 = axes('Position', [0.30 0.25 0.05 .18]);

    rho = []; p = [];
    for indSub = 1:32
        xData = Figure8.e_abn_singleTrial(1,indSub,:,:);
        yData = Figure8.eDiffRel_singleTrial(1,indSub,:,:);
        [rho(indSub), p(indSub)] = corr(xData(:),yData(:), 'type', 'Pearson', 'rows', 'complete');
        z(indSub) = .5*log((1+rho(indSub))/(1-rho(indSub))); % transform to Fisher's z
    end

    hi = histogram(z,[-1:.2:1]); % Note that across individuals there are only positive associations within subject!
    hi.FaceColor = 'k';
    title({'Within-subject';'correlations'});
    xlabel('Fisher`s z'); ylabel('# Subjects');
    ylim([0 15])

    set(handaxes1, 'Box', 'off')
    % Adjust XY label font
    set(get(handaxes1, 'XLabel'), 'FontName', 'Times')
    set(get(handaxes1, 'YLabel'), 'FontName', 'Times')
    set(findall(handaxes1,'-property','FontSize'),'FontSize',15)

    [h, p_2ndLevel_eSNR,ci,stats] = ttest(z)
    nanmean(z)

    if p_2ndLevel_eSNR==0
        p_2ndLevel_eSNR = 1^-14;
    end
    
handaxes2 = axes('Position', [0.59 0.25 0.05 .18]);

    rho = []; p = [];
    for indSub = 1:32
        xData = Figure8.e_abn_singleTrial(1,indSub,:,:);
        yData = Figure8.e_BG_singleTrial(1,indSub,:,:);
        [rho(indSub), p(indSub)] = corr(xData(:),yData(:), 'type', 'Pearson', 'rows', 'complete');
        z(indSub) = .5*log((1+rho(indSub))/(1-rho(indSub))); % transform to Fisher's z
    end

    hi = histogram(z,[-1:.2:1]); % Note that across individuals there are only positive associations within subject!
    hi.FaceColor = 'k';
    title({'Within-subject';'correlations'});
    xlabel('Fisher`s z'); ylabel('# Subjects');

    ylim([0 15])
    set(handaxes2, 'Box', 'off')
    % Adjust XY label font
    set(get(handaxes2, 'XLabel'), 'FontName', 'Times')
    set(get(handaxes2, 'YLabel'), 'FontName', 'Times')
    set(findall(handaxes2,'-property','FontSize'),'FontSize',15)
    
    [h, p_2ndLevel_BG] = ttest(z)
    nanmean(z)
    
pn.plotFolder = '/Volumes/EEG/BOSC/BOSC_SternRest/X_documentation/B_2018_Manuscript/Figures_resubmit1/F1_Figures/';
figureName = 'F7AB';

saveas(h, [pn.plotFolder, figureName], 'fig');
saveas(h, [pn.plotFolder, figureName], 'epsc');
saveas(h, [pn.plotFolder, figureName], 'png');
