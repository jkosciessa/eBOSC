% Calculate amplitude-abundance link across load levels

% based originally on /Volumes/EEG/BOSC_Sternberg/scripts/B_BOSC/C_JQKplotAmpAbn_170703.mlx

% 180920 | use updated X structure with correct alpha range

% %% set paths
% 
% pn.root = '/Volumes/EEG/BOSC_Sternberg/B_analyses/B_BOSC/B_eBOSC_170703_full_wl6/';
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
% save('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/F8.mat', 'Figure8')

%% load Figure data

load('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/F8.mat', 'Figure8')

%% Plot with inlays
% see https://de.mathworks.com/matlabcentral/fileexchange/35245-matlab-plot-gallery-plot-in-plot?focused=6792989&tab=example
    
h = figure('DefaultAxesFontSize',12, 'units','normalized','position',[.1 .1 .7 .8]);
handaxes1 = axes('Position', [0.12 0.55 0.8 0.4]);

    for indSub = 1:32
        xData = Figure8.e_abn_singleTrial(1,indSub,:,:);
        yData = Figure8.eDiffRel_singleTrial(1,indSub,:,:);
        scatter(xData(:),yData(:),[],Figure8.colormWithinSub(indSub,:),'filled');
        hold on;
    end
    title('Single-trial estimates, color-coded by subject', 'FontSize', 20);
    xlabel('Abundance'); ylabel ('Rhythmic SNR [a.u.]');
    xlim([0 1]);

    set(handaxes1, 'Box', 'off')
    % Adjust XY label font
    handxlabel1 = get(gca, 'XLabel');
    set(handxlabel1, 'FontSize', 16, 'FontWeight', 'bold')
    handylabel1 = get(gca, 'ylabel');
    set(handylabel1, 'FontSize', 16, 'FontWeight', 'bold')

% add histogram for individual amplitude-abundance correlations

handaxes2 = axes('Position', [0.78 0.60 0.11 0.13]);

    rho = []; p = [];
    for indSub = 1:32
        xData = Figure8.e_abn_singleTrial(1,indSub,:,:);
        yData = Figure8.eDiffRel_singleTrial(1,indSub,:,:);
        [rho(indSub), p(indSub)] = corr(xData(:),yData(:), 'type', 'Pearson', 'rows', 'complete');
        z(indSub) = .5*log((1+rho(indSub))/(1-rho(indSub))); % transform to Fisher's z
    end

    % compute the conditional association
    [r, p] = corrcoef(z, squeeze(nanmean(nanmean(Figure8.eDiffRel_singleTrial,4),3)))
    
% r =
% 
%     1.0000   -0.0724
%    -0.0724    1.0000
% 
% 
% p =
% 
%     1.0000    0.6936
%     0.6936    1.0000

    
    hi = histogram(z,[-1:.2:1]); % Note that across individuals there are only positive associations within subject!
    hi.FaceColor = 'k';
    title('Within-subject correlations');
    xlabel('Fisher`s z'); ylabel('Subject amount');

    set(handaxes2, 'Box', 'off')
    % Adjust XY label font
    set(get(handaxes2, 'XLabel'), 'FontName', 'Times')
    set(get(handaxes2, 'YLabel'), 'FontName', 'Times')
    
%% Background correlations

%h = figure('DefaultAxesFontSize',12, 'units','normalized','position',[.1 .1 .7 .4]);
handaxes1 = axes('Position', [0.12 0.06 0.8 0.4]);

    for indSub = 1:32
        xData = Figure8.e_abn_singleTrial(1,indSub,:,:);
        yData = Figure8.e_BG_singleTrial(1,indSub,:,:);
        scatter(xData(:),yData(:),[],Figure8.colormWithinSub(indSub,:),'filled');
        hold on;
    end
    %title('Single-trial estimates, color-coded by subject', 'FontSize', 20);
    xlabel('Abundance'); ylabel ('Background amplitude [a.u.]');
    xlim([0 1]);

    set(handaxes1, 'Box', 'off')
    % Adjust XY label font
    handxlabel1 = get(gca, 'XLabel');
    set(handxlabel1, 'FontSize', 16, 'FontWeight', 'bold')
    handylabel1 = get(gca, 'ylabel');
    set(handylabel1, 'FontSize', 16, 'FontWeight', 'bold')

% add histogram for individual amplitude-abundance correlations

handaxes2 = axes('Position', [0.78 0.11 0.11 .13]);

    rho = []; p = [];
    for indSub = 1:32
        xData = Figure8.e_abn_singleTrial(1,indSub,:,:);
        yData = Figure8.e_BG_singleTrial(1,indSub,:,:);
        [rho(indSub), p(indSub)] = corr(xData(:),yData(:), 'type', 'Pearson', 'rows', 'complete');
        z(indSub) = .5*log((1+rho(indSub))/(1-rho(indSub))); % transform to Fisher's z
    end

    hi = histogram(z,[-1:.2:1]); % Note that across individuals there are only positive associations within subject!
    hi.FaceColor = 'k';
    title('Within-subject correlations');
    xlabel('Fisher`s z'); ylabel('Subject amount');

    set(handaxes2, 'Box', 'off')
    % Adjust XY label font
    set(get(handaxes2, 'XLabel'), 'FontName', 'Times')
    set(get(handaxes2, 'YLabel'), 'FontName', 'Times')
    
pn.plotFolder = '/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F_Figures/';
figureName = 'F8';

saveas(h, [pn.plotFolder, figureName], 'fig');
saveas(h, [pn.plotFolder, figureName], 'epsc');
saveas(h, [pn.plotFolder, figureName], 'png');