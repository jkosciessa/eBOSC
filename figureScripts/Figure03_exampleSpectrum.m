% Plot individual rhythm-conditional spectra
% adapted from v3 containing the computations

%% load output

clear all; clc;

pn.root = '/Users/kosciessa/Desktop/mntTardis/Stern_WIP/B_eBOSC_180917_retention/';
pn.spectraOut = [pn.root, 'X_spectra/B_data/'];
pn.plotFolder = [pn.root, 'X_spectra/C_figures/']; if ~exist(pn.plotFolder); mkdir(pn.plotFolder); end

load([pn.spectraOut, 'A_meanSpectra.mat'], 'results');

freq = 2.^[0:.125:6];

%% Plot example sepctrum with rhythm-conditional variant

channels = 44:60;
freq = 2.^[0:.125:6];
h = figure('units','normalized','position',[.1 .1 .25 .3]);
hold on;
plot(squeeze(nanmedian(nanmedian(results{3}.L6o(4,channels,:),2),1)), 'Color', [0 0 0],'LineWidth', 2)
plot(squeeze(nanmedian(nanmedian(results{3}.L6a(4,channels,:),2),1)), 'Color', [1 0 0],'LineWidth', 2)
xticks([1:5:49]); xticklabels(round(freq(xticks),1)); xlabel('Frequency')
ylabel('Amplitude')
xlim([1 49])
set(findall(gcf,'-property','FontSize'),'FontSize',18)

pn.plotFolder = '/Users/kosciessa/Desktop/';
figureName = 'A_exampleSpectrum';

%saveas(h, [pn.plotFolder, figureName], 'fig');
saveas(h, [pn.plotFolder, figureName], 'epsc');
%saveas(h, [pn.plotFolder, figureName], 'png');