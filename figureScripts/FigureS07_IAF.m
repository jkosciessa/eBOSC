% original: /Volumes/EEG/BOSC_Sternberg/scripts/B_BOSC/K_IAFVariability_170830.m

% check stabilty of single-trial abundance, amplitude

% 170830: initially wrong referencing of sessions (e.g. CuttingEEG); now corrected

% clear all; clc;
% 
% pn.root = '/Volumes/EEG/BOSC_Sternberg/B_analyses/B_BOSC/B_eBOSC_170703_full_wl6/';
% pn.boscData = [pn.root, 'B_extractIndices/B_data/'];
% pn.simData = '/Volumes/EEG/BOSC_Simulation/17_JQK_RED/D_data/RED_161212/';
% pn.oIAFdata = [pn.root, 'K_IAFVariability/B_data/'];
% 
% info.channels = 44:60; % restrict to parieto-occipital channels
% cond = {'L2'; 'L4'; 'L6'};
% 
% %% load data
% 
% load([pn.boscData, 'X_8to15_170703_v3.mat']);
% load([pn.simData, 'IAFResults.mat'], 'SimData');
% Figure10.SimData = SimData;
% load([pn.oIAFdata, 'A_IAF_ByTrial_170703_8_15.mat'], 'oIAF_byTrial');
% 
% %% extract data
% 
% for indSession = 1:3
%     Figure10.mean_abn{indSession} = []; Figure10.std_iaf_A{indSession} = [];
%     for indCond = 1:3
%         abundance = X{indSession,1}.([cond{indCond},'e_abn']); % sub*chan*trial
%         chanAvg_abn = squeeze(nanmean(abundance(:,info.channels,:),2));
%         Figure10.mean_abn{indSession} = [Figure10.mean_abn{indSession}; nanmean(chanAvg_abn,2)];
% 
%         iaf = X{indSession,1}.([cond{indCond},'e_iaf']); % sub*chan*trial
%         chanAvg_iaf = squeeze(nanmean(iaf(:,info.channels,:),2));
%         Figure10.std_iaf_A{indSession} = [Figure10.std_iaf_A{indSession}; nanstd(chanAvg_iaf,[],2)];
%     end
% end
% 
% for indSession = 1:3
%     iaf = NaN(32,60,70);
%     for indSub = 1:32
%         amountTrials = size(oIAF_byTrial{indSession,indSub},2);
%         iaf(indSub,:,1:amountTrials) = oIAF_byTrial{indSession,indSub}; % create sub*channel*trial matrix
%     end
%     Figure10.mean_snr{indSession} = []; Figure10.std_iaf_B{indSession} = [];
%     for indCond = 1:3
%         snr = repmat(X{indSession,1}.([cond{indCond},'o_amp'])./nanmean(X{indSession,1}.([cond{indCond},'e_fitBG']),3),1,1,70);
%         chanAvg_snr = squeeze(nanmean(snr(:,info.channels,:),2));
%         Figure10.mean_snr{indSession} = [Figure10.mean_snr{indSession}; nanmean(chanAvg_snr,2)];
% 
%         chanAvg_iaf = squeeze(nanmean(iaf(:,info.channels,:),2));
%         Figure10.std_iaf_B{indSession} = [Figure10.std_iaf_B{indSession}; nanstd(chanAvg_iaf,[],2)];
%     end
% end
% 
% %% save Figure data
% 
% save('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/F10.mat', 'Figure10')

%% load Figure data

load('/Users/kosciessa/Desktop/eBOSC/figureData/S7.mat', 'Figure10')

%% Figure on IAF variability

h = figure('units','normalized','position',[.1 .1 .7 .3]);
subplot(1,2,1);
    for indSession = 1:3
        line{indSession} = scatter(Figure10.mean_abn{indSession}, ...
            Figure10.std_iaf_A{indSession}, 'filled'); hold on;
        r{indSession} = corrcoef(Figure10.mean_abn{indSession}, ...
            Figure10.std_iaf_A{indSession}); % -.8 correlation between abundance mean across trials and IAF variability
    end
    xlabel('Abundance'); ylabel('rhythmic IAF variability [Hz]');
    xlim([0 1]);
    title('Rhythmic IAF variability across trials')
    set(findall(gcf,'-property','FontSize'),'FontSize',18)

    % overlay simulation results 

    hold on;
    plot(Figure10.SimData.abn_AcrossDurs, Figure10.SimData.stdIAF_AcrossDurs, ...
        'k', 'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 2);

    legend([line{1}, line{2}, line{3}], {['Session 1: r = ', sprintf('%.2f', r{1}(2)), ', p < .001']; ...
        ['Session 7: r = ', sprintf('%.2f', r{2}(2)), ', p < .001']; ...
        ['Session 8: r = ', sprintf('%.2f', r{3}(2)), ', p < .001']}); legend('boxoff');
    
subplot(1,2,2);

    % load IAF trial-by-trial estimates from the whole-trial signal (i.e. overall vs. constrained to rhythmic timepoints)

    for indSession = 1:3
        line{indSession} = scatter(Figure10.mean_snr{indSession}, Figure10.std_iaf_B{indSession}, 'filled'); hold on;
        availableTrials = find(~isnan(Figure10.mean_snr{indSession}));
        [r{indSession}, p{indSession}] = corrcoef(Figure10.mean_snr{indSession}(availableTrials), ...
            Figure10.std_iaf_B{indSession}(availableTrials));
    end
    xlabel('Overall SNR'); ylabel('whole-trial IAF variability [Hz]');
    xlim([0 10]);
    title('Whole-trial IAF variability across trials')

    hold on;
    plot(Figure10.SimData.meanSNR_overall_AcrossDurs,...
        Figure10.SimData.stdIAF_overall_AcrossDurs, ...
        'k', 'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 1);

    legend([line{1}, line{2}, line{3}],{['Session 1: r = ', sprintf('%.2f', r{1}(2)), ', p < .001']; ...
        ['Session 7: r = ', sprintf('%.2f', r{2}(2)), ', p < .001']; ...
        ['Session 8: r = ', sprintf('%.2f', r{3}(2)), ', p < .001']}); legend('boxoff');
    
    set(findall(gcf,'-property','FontSize'),'FontSize',18)

%% save Figure
    
pn.plotFolder = ['/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F_Figures/'];
figureName = 'F10';

saveas(h, [pn.plotFolder, figureName], 'fig');
saveas(h, [pn.plotFolder, figureName], 'epsc');
saveas(h, [pn.plotFolder, figureName], 'png');
