%% Plot rhythm-conditional spectra
% 
% pn.root = '/Volumes/EEG/BOSC_Sternberg/B_analyses/B_BOSC/B_eBOSC_180917_retention/';
% pn.spectraOut = [pn.root, 'X_spectra/B_data/'];
% 
% load([pn.spectraOut, 'A_meanSpectra.mat'], 'results');
% 
% % average across sessions and loads
% 
% Figure12A.o_across = cat(4, results{1}.L2o, results{1}.L4o, results{1}.L6o, ...
%     results{2}.L2o, results{2}.L4o, results{2}.L6o, ...
%     results{3}.L2o, results{3}.L4o, results{3}.L6o);
% Figure12A.o_across = squeeze(mean(Figure12A.o_across,4));
% 
% Figure12A.nt_across = cat(4, results{1}.L2nt, results{1}.L4nt, results{1}.L6nt, ...
%     results{2}.L2nt, results{2}.L4nt, results{2}.L6nt, ...
%     results{3}.L2nt, results{3}.L4nt, results{3}.L6nt);
% Figure12A.nt_across = squeeze(mean(Figure12A.nt_across,4));
% 
% Figure12A.na_across = cat(4, results{1}.L2na, results{1}.L4na, results{1}.L6na, ...
%     results{2}.L2na, results{2}.L4na, results{2}.L6na, ...
%     results{3}.L2na, results{3}.L4na, results{3}.L6na);
% Figure12A.na_across = squeeze(mean(Figure12A.na_across,4));
% 
% Figure12A.nb_across = cat(4, results{1}.L2nb, results{1}.L4nb, results{1}.L6nb, ...
%     results{2}.L2nb, results{2}.L4nb, results{2}.L6nb, ...
%     results{3}.L2nb, results{3}.L4nb, results{3}.L6nb);
% Figure12A.nb_across = squeeze(mean(Figure12A.nb_across,4));
% 
% Figure12A.ng_across = cat(4, results{1}.L2ng, results{1}.L4ng, results{1}.L6ng, ...
%     results{2}.L2ng, results{2}.L4ng, results{2}.L6ng, ...
%     results{3}.L2ng, results{3}.L4ng, results{3}.L6ng);
% Figure12A.ng_across = squeeze(mean(Figure12A.ng_across,4));
% 
% Figure12A.t_across = cat(4, results{1}.L2t, results{1}.L4nt, results{1}.L6t, ...
%     results{2}.L2t, results{2}.L4t, results{2}.L6t, ...
%     results{3}.L2t, results{3}.L4t, results{3}.L6t);
% Figure12A.t_across = squeeze(nanmean(Figure12A.t_across,4));
% 
% Figure12A.a_across = cat(4, results{1}.L2a, results{1}.L4a, results{1}.L6a, ...
%     results{2}.L2a, results{2}.L4a, results{2}.L6a, ...
%     results{3}.L2a, results{3}.L4a, results{3}.L6a);
% Figure12A.a_across = squeeze(nanmean(Figure12A.a_across,4));
% 
% Figure12A.b_across = cat(4, results{1}.L2b, results{1}.L4b, results{1}.L6b, ...
%     results{2}.L2b, results{2}.L4b, results{2}.L6b, ...
%     results{3}.L2b, results{3}.L4b, results{3}.L6b);
% Figure12A.b_across = squeeze(nanmean(Figure12A.b_across,4));
% 
% Figure12A.g_across = cat(4, results{1}.L2g, results{1}.L4g, results{1}.L6g, ...
%     results{2}.L2g, results{2}.L4g, results{2}.L6g, ...
%     results{3}.L2g, results{3}.L4g, results{3}.L6g);
% Figure12A.g_across = squeeze(nanmean(Figure12A.g_across,4));
% 
% Figure12A.freq = 2.^[0:.125:6];

% %% save Figure data
% 
% save('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/F12A.mat', 'Figure12A')

%% load Figure data

load('/Users/kosciessa/Desktop/eBOSC/figureData/F9A.mat', 'Figure12A')

pn.shadedError = ['/Volumes/EEG/BOSC/BOSC_Sternberg/T_tools/shadedErrorBar-7002ebc']; addpath(pn.shadedError);

%% Plot rhythm-conditional deviance spectra

h = figure('units','normalized','position',[.1 .1 .5 .5]);
hold on; 
curMeas = nanmean(Figure12A.t_across-Figure12A.nt_across,2);
    standError = nanstd(curMeas,1)./sqrt(size(curMeas,1));
    l1 = shadedErrorBar([],nanmean(curMeas,1),standError, 'lineprops', {'linewidth', 4}, 'patchSaturation', .05);
curMeas = nanmean(Figure12A.a_across-Figure12A.na_across,2);
    standError = nanstd(curMeas,1)./sqrt(size(curMeas,1));
    l2 = shadedErrorBar([],nanmean(curMeas,1),standError, 'lineprops', {'linewidth', 4}, 'patchSaturation', .05);
curMeas = nanmean(Figure12A.b_across-Figure12A.nb_across,2);
    standError = nanstd(curMeas,1)./sqrt(size(curMeas,1));
    l3 = shadedErrorBar([],nanmean(curMeas,1),standError, 'lineprops', {'linewidth', 4}, 'patchSaturation', .05);
curMeas = nanmean(Figure12A.g_across-Figure12A.ng_across,2);
    standError = nanstd(curMeas,1)./sqrt(size(curMeas,1));
    l4 = shadedErrorBar([],nanmean(curMeas,1),standError, 'lineprops', {'linewidth', 4}, 'patchSaturation', .05);
curMeas = nanmean(Figure12A.o_across,2);
    standError = nanstd(curMeas,1)./sqrt(size(curMeas,1));
    l5 = shadedErrorBar([],nanmean(curMeas,1),standError, 'lineprops', {'Color', [0 0 0], 'linewidth', 8}, 'patchSaturation', .05);
ax = gca; ax.XTick = [1:4:numel(Figure12A.freq)]; ax.XTickLabels = round(Figure12A.freq(ax.XTick),1); xlim([1 numel(Figure12A.freq)])
xlabel('Frequency [Hz]'); ylabel('Amplitude [a.u.]');
legend([l5.mainLine, l1.mainLine,l2.mainLine,l3.mainLine,l4.mainLine], ...
    {'Overall spectrum'; 'Theta-conditional spectrum'; 'Alpha-conditional spectrum';...
    'Beta-conditional spectrum'; 'Gamma-conditional spectrum'});
legend('boxoff');
title([{'Temporal knowledge about rhythmic episodes';'can be used to obtain rhythm-conditional spectra'}]);
set(findall(gcf,'-property','FontSize'),'FontSize',30)

pn.plotFolder = '/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F_Figures/';
figureName = 'F12A';

saveas(h, [pn.plotFolder, figureName], 'fig');
saveas(h, [pn.plotFolder, figureName], 'epsc');
saveas(h, [pn.plotFolder, figureName], 'png');
