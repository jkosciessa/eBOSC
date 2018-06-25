% original: /Volumes/EEG/BOSC_Simulation/17_JQK_RED/A_scripts/HX_FigureAmpAbn_JQK_170630_wOverall_v4.m

%% 
% This script plots the results from the siumulation analysis and pits the 
% standard BOSC procedure against two variants of extended BOSC. 
% 
% THG denotes Thomas' post-processing for wavelet convolution; FWHM denotes 
% the full-width at half maximum approach. Otherwise, the extended BOSC variants 
% are the same, making the detected matrix sparser by finding consecutive rhythmic 
% 'episodes' and then conducting the correction for wavelet smearing. Beforehand, 
% the background has been fit with undersampling in the alpha range and robust 
% regression for the extended BOSC variant. Standard BOSC uses linear regression 
% with the full spectrum.
% 
% Hits and FAs are calculated from the time series indices of detected alpha 
% episodes (in comparison to the indices of simulated alpha). As no episodes are 
% available for standard BOSC, Hits and FAs have been computed considering every 
% sample point with a detected point in the alpha range (8-12 Hz) from the full 
% detected matrix. Abundance calculated from this vector (i.e. pepisode of all 
% detected point in alpha range) is denoted as 'abundance_PepMeanAlpha'. Alternatively, 
% the average of pepisode values is denoted as 'abundance_meanPep'.
% 
% 
% 
% Simulation parameters:
% 
% - 9.5 Hz simulated with varying amplitudes and cycles (see plots for values) 
% centered in 20s time interval, added on 1/f background
% 
% - sampling rate = 250 Hz; 500 samples removed at each edge for artifacts, 
% 250 sampled removed from BOSC detection matrix
% 
% - effective 14 s interval
% 
% - 200 iterations (save computing time for now) with BG calculated only 
% from the first 100 trials

pn.data = '/Volumes/EEG/BOSC_Simulation/17_JQK_RED/D_data/D_simulation/';
pn.plotFolder = '/Volumes/EEG/BOSC_Simulation/17_JQK_RED/B_figures/';

% method x amplitude x condition x repetitions [for extended BOSC]
% condition: 1-pure sine; 2- sine with BG

addpath(genpath('/Volumes/EEG/BOSC_SternRest/P&G/X_Other/DrosteEffect-BrewerMap-7fdcfcf/'))
colorm = brewermap(4,'RdYlBu'); % use external function to create colormap

%% Load simulation results

simulationSetup = 4;

load([pn.data, 'AmpAbnSimulation_v4_170630.mat']);

SNR = squeeze(SignalDetection.Amp(:,:,simulationSetup,:))./squeeze(SignalDetection.fitBG(:,:,simulationSetup,:));
empiricalSNR = round(nanmean(SNR,2),1);
for a = 1:numel(amplitude)
    ampSNR{a} = [num2str(amplitude(a)), ' (', num2str(empiricalSNR(a)),') '];
end

oAmp = squeeze(SignalDetection.oAmp(:,:,simulationSetup,:));
oBG = squeeze(SignalDetection.oBGestimate(:,:,simulationSetup,:));

SNRo = squeeze(SignalDetection.oAmp(:,:,simulationSetup,:))./squeeze(SignalDetection.oBGestimate(:,:,simulationSetup,:));
empiricalSNRo = round(nanmean(SNRo(:,:),2),1);
for a = 1:numel(amplitude)
    ampSNR{a} = [num2str(amplitude(a)), ' (', num2str(empiricalSNRo(a)),') '];
end

%% 

h = figure('units','normalized','position',[.1 .1 .7 .4]);
subplot(1,2,1);
scatter(mean(SignalDetection.abn(:,:,4,:),4), mean(SignalDetection.Amp(:,:,2,:),4),60,'MarkerEdgeColor', [1,1,1], 'MarkerFaceColor', [0,0,0]); hold on; % simulated amplitudes
scatter(mean(SignalDetection.abn(:,:,4,:),4), mean(SignalDetection.Amp(:,:,4,:),4),60, 'MarkerEdgeColor', [1,1,1], 'MarkerFaceColor', colorm(1,:)); % detected rhythmic amplitudes
scatter(mean(SignalDetection.abn(:,:,4,:),4), mean(SignalDetection.Amp(:,:,4,:),4)-mean(SignalDetection.fitBG(:,:,4,:),4),60, 'MarkerEdgeColor', [1,1,1], 'MarkerFaceColor', colorm(2,:)); % detected rhythmic amplitude-BG
% add overall amplitude
scatter(mean(SignalDetection.abn(:,:,4,:),4), mean(SignalDetection.oAmp(:,:,4,:),4),60,'MarkerEdgeColor', [1,1,1], 'MarkerFaceColor', colorm(4,:)); hold on; % overall amplitudes
ylabel('Empirical Amplitude [a.u.]'); xlabel('Empirical Abundance');
%title('Simulated vs. empirical amplitudes');
legend({'simulated data'; 'rhythmic estimate'; ['rhythmic estimate' char(10) 'excl. background']; 'overall estimate'}, 'Location', 'northwest')
legend('boxoff')
legend('AutoUpdate','off')
title('Simulated Abundance: 1');
line([1 1], [0, 150], 'Color','k', 'LineWidth', 2) % 32 cycles in 14 seconds

% fit to lower simulated abundance

subplot(1,2,2);
scatter(mean(SignalDetection.abn(:,:,3,:),4), mean(SignalDetection.Amp(:,:,1,:),4),60,'MarkerEdgeColor', [1,1,1], 'MarkerFaceColor', [0,0,0]); hold on; % simulated amplitudes
scatter(mean(SignalDetection.abn(:,:,3,:),4), mean(SignalDetection.Amp(:,:,3,:),4),60, 'MarkerEdgeColor', [1,1,1], 'MarkerFaceColor', colorm(1,:)); % detected rhythmic amplitudes
scatter(mean(SignalDetection.abn(:,:,3,:),4), mean(SignalDetection.Amp(:,:,3,:),4)-mean(SignalDetection.fitBG(:,:,3,:),4),60, 'MarkerEdgeColor', [1,1,1], 'MarkerFaceColor', colorm(2,:)); % detected rhythmic amplitude-BG
scatter(mean(SignalDetection.abn(:,:,3,:),4), mean(SignalDetection.oAmp(:,:,3,:),4),60,'MarkerEdgeColor', [1,1,1], 'MarkerFaceColor', colorm(4,:)); hold on; % overall amplitudes
ylabel('Empirical Amplitude [a.u.]'); xlabel('Empirical Abundance');
title('Simulated Abundance: .25');

set(findall(gcf,'-property','FontSize'),'FontSize',16)

% add SNR as text
x = mean(SignalDetection.abn(:,:,1,:),4);
y = mean(SignalDetection.Amp(:,:,1,:),4);
z = ampSNR;
for K = 1:size(mean(SignalDetection.abn(:,:,3,:),4),2)
    text(x(K)+.02,y(K),z{K}, 'FontSize', 13)
end

line([(35/10)/14 (35/10)/14], [0, 150], 'Color','k', 'LineWidth', 2) % 32 cycles in 14 seconds

saveas(h, [pn.plotFolder, datestr(now, 'yymmdd'), '_SimulatedAmplitude'], 'fig');
saveas(h, [pn.plotFolder, datestr(now, 'yymmdd'), '_SimulatedAmplitude'], 'epsc');

%% save amp-abn relationship for comparison with empirical resting state

figure;
plot(mean(SignalDetection.abn(:,:,4,:),4), empiricalSNR);

simData.Abn = mean(SignalDetection.abn(:,:,4,:),4);
simData.SNR = empiricalSNR;

save([pn.data, '_SimDataforRestingStateComparison_v4.mat'], 'simData');

figure;
plot(mean(SignalDetection.abn(:,:,4,:),4), empiricalSNRo);

simData.Abn = mean(SignalDetection.abn(:,:,4,:),4);
simData.SNR = empiricalSNRo;

save([pn.data, '_SimDataforRestingStateComparison_overall_v4.mat'], 'simData');