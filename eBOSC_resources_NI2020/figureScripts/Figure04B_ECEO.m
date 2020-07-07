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

load('/Users/kosciessa/Desktop/eBOSC/figureData/F4B.mat', 'Figure5B')

%% Figure: Eye closure state change

addpath('/Volumes/EEG/BOSC/BOSC_SternRest/X_documentation/B_2018_Manuscript/F4_FigureTools/')

addpath('/Volumes/EEG/BOSC/BOSC_SternRest/B_analyses/A_eBOSC/A_SternBRest_170703/C_interrelations/T_tools/brewermap/')
Figure5B.colorm = brewermap(4,'RdYlBu'); % use external function to create colormap

h = figure('units','normalized','position',[.1 .1 .25 .7]);
subplot(2,1,1); hold on;
    x = Figure5B.EC1e_abn; y = log10(Figure5B.EC1e_amp); scatter(x, y, 75,'filled', 'MarkerFaceColor', Figure5B.colorm(1,:)); 
    y1_ls = polyval(polyfit(x,y,1),x); y1_ls = plot(x, y1_ls, 'Color', Figure5B.colorm(1,:), 'LineWidth', 4);
    [rl1,pl1] = corrcoef(x,y); pval1 = []; pval1 = convertPtoExponential(pl1(2));
    x = Figure5B.EC2e_abn; y = log10(Figure5B.EC2e_amp); scatter(x, y, 75,'filled', 'MarkerFaceColor', Figure5B.colorm(2,:)); 
    y2_ls = polyval(polyfit(x,y,1),x); y2_ls = plot(x, y2_ls, 'Color', Figure5B.colorm(2,:), 'LineWidth', 4);
    [rl2,pl2] = corrcoef(x,y); pval2 = []; pval2 = convertPtoExponential(pl2(2));
    x = Figure5B.EO1e_abn; y = log10(Figure5B.EO1e_amp); scatter(x, y, 75,'filled', 'MarkerFaceColor', Figure5B.colorm(3,:)); 
    y3_ls = polyval(polyfit(x,y,1),x); y3_ls = plot(x, y3_ls, 'Color', Figure5B.colorm(3,:), 'LineWidth', 4);
    [rl3,pl3] = corrcoef(x,y); pval3 = []; pval3 = convertPtoExponential(pl3(2));
    x = Figure5B.EO2e_abn; y = log10(Figure5B.EO2e_amp); scatter(x, y, 75,'filled', 'MarkerFaceColor', Figure5B.colorm(4,:)); 
    y4_ls = polyval(polyfit(x,y,1),x); y4_ls = plot(x, y4_ls, 'Color', Figure5B.colorm(4,:), 'LineWidth', 4);
    [rl4,pl4] = corrcoef(x,y); pval4 = []; pval4 = convertPtoExponential(pl4(2));
    leg1 = legend([y1_ls, y2_ls, y3_ls, y4_ls], {['EC: continuous: r=',num2str(round(rl1(2),2)), ', p=' pval1{1}];...
        ['EC: interleaved: r=',num2str(round(rl2(2),2)), ', p=' pval2{1}];...
        ['EO: continuous: r=',num2str(round(rl3(2),2)), ', p=' pval3{1}];...
        ['EO: interleaved: r=',num2str(round(rl4(2),2)), ', p=' pval4{1}]}, 'location', 'SouthEast'); legend('boxoff');
    xlabel('Abundance'); ylabel('Rhythmic Amplitude (log10)');
    ylim([1.8 3.4])
    title({'Amplitude and abundance are related within condition'; ''})
subplot(2,1,2); hold on;
    x = Figure5B.EC1e_abn-Figure5B.EO1e_abn; y = log10(Figure5B.EC1e_amp)-log10(Figure5B.EO1e_amp); scatter(x, y, 75,'filled', 'MarkerFaceColor', Figure5B.colorm(1,:)); 
    y1_ls = polyval(polyfit(x,y,1),x); y1_ls = plot(x, y1_ls, 'Color', Figure5B.colorm(1,:), 'LineWidth', 4);
    [rl1,pl1] = corrcoef(x,y); pval1 = []; pval1 = convertPtoExponential(pl1(2));
    x = Figure5B.EC2e_abn-Figure5B.EO2e_abn; y = log10(Figure5B.EC2e_amp)-log10(Figure5B.EO2e_amp); scatter(x, y, 75,'filled', 'MarkerFaceColor', Figure5B.colorm(4,:)); 
    y2_ls = polyval(polyfit(x,y,1),x); y2_ls = plot(x, y2_ls, 'Color', Figure5B.colorm(4,:), 'LineWidth', 4);
    [rl2,pl2] = corrcoef(x,y); pval2 = []; pval2 = convertPtoExponential(pl2(2));
    xlabel('Difference in Abundance'); ylabel('Difference in Rhythmic Amplitude (log10)');
    leg2 = legend([y1_ls, y2_ls], {['EC-EO: continuous: r=',num2str(round(rl1(2),2)), ', p=' pval1{1}];...
        ['EC-EO: interleaved: r=',num2str(round(rl2(2),2)), ', p=' pval2{1}]}, 'location', 'SouthEast'); legend('boxoff');
        title({'Berger effect of amplitude and abundance are collinear'; ''})
        ylim([-.2 .8]); xlim([-.2 1])
    set(findall(gcf,'-property','FontSize'),'FontSize',21)

    set(leg1, 'FontSize', 18)
    set(leg2, 'FontSize', 18)

    
% h = figure('units','normalized','position',[.1 .1 .3 .7]);
% title('Parameter change upon eye closure (individual subjects)')
% for j = 1:32 % EO1-EC1
%     hold on; h1 = arrow([Figure5B.EO1e_abn(j) Figure5B.EO1e_amp(j)],[Figure5B.EC1e_abn(j) Figure5B.EC1e_amp(j)], 5, 'Color',Figure5B.colorm(1,:));
% end; clear j
% hold on; arrow([mean(Figure5B.EO1e_abn) mean(Figure5B.EO1e_amp)],[mean(Figure5B.EC1e_abn) mean(Figure5B.EC1e_amp)], 5, 'Color',[.5 0 0], 'LineWidth', 3);
% for j = 1:32 % EO2-EC2
%     hold on; h2 = arrow([Figure5B.EO2e_abn(j) Figure5B.EO2e_amp(j)],[Figure5B.EC2e_abn(j) Figure5B.EC2e_amp(j)], 5, 'Color',Figure5B.colorm(4,:));
% end; clear j
% hold on; arrow([mean(Figure5B.EO2e_abn) mean(Figure5B.EO2e_amp)],[mean(Figure5B.EC2e_abn) mean(Figure5B.EC2e_amp)], 5, 'Color',[0 0 .5], 'LineWidth', 3);
% xlabel('Abundance'); ylabel('Rhythmic Amplitude excl. BG [?V]');
% legend([h1, h2], {'continuous', 'interleaved'}, 'location', 'northwest');
% legend('boxoff')
% set(findall(gcf,'-property','FontSize'),'FontSize',15)

pn.plotFolder = '/Volumes/EEG/BOSC/BOSC_SternRest/X_documentation/B_2018_Manuscript/Figures_resubmit1/F1_Figures/';

saveas(h, [pn.plotFolder, 'F4B'], 'fig');
saveas(h, [pn.plotFolder, 'F4B'], 'epsc');
