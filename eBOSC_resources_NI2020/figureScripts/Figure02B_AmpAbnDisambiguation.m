%% Plot disambiguation of rhythmic power and duration (here exemplified by alpha)

% %% create Figure data
% 
% Sine = load(['/Volumes/EEG/BOSC_Simulation/17_JQK_RED/A_scripts/N_pureSine_AmpAbn/B_data/AmpAbnSimulation_pureSineOnly.mat']);
% ExtendedB = load(['/Volumes/EEG/BOSC_Simulation/17_JQK_RED/D_data/D_simulation/REDSimulation_171023_v10B.mat']);
% 
% Figure3B.Sine_oAmp = squeeze(Sine.SignalDetection.oAmp(:,2:end,1:14));
% Figure3B.Sine_Amp = squeeze(Sine.SignalDetection.Amp(:,2:end,1:14));
% Figure3B.Sine_Abn = squeeze(Sine.SignalDetection.abn(:,2:end,1:14));
% Figure3B.Sine_Ratio = squeeze(Sine.SignalDetection.oAmp(:,2:end,1:14)./...
%     Sine.SignalDetection.Amp(:,2:end,1:14));
% 
% Figure3B.Signal_oAmp = squeeze(nanmean(ExtendedB.SignalDetection.oAmp(2:end,2:end,:),3)); % average across trials
% Figure3B.Signal_Amp = squeeze(nanmean(ExtendedB.SignalDetection.Amp(2:end,2:end,:),3));
% Figure3B.Signal_Abn = squeeze(nanmean(ExtendedB.SignalDetection.abn(2:end,2:end,:),3));
% Figure3B.Signal_Ratio = squeeze(nanmean(ExtendedB.SignalDetection.oAmp(2:end,2:end,:)./...
%     ExtendedB.SignalDetection.Amp(2:end,2:end,:),3));
% 
% for c = 1:numel(Sine.cfg.simParams.cycles)
%     Abundance(1,c) = Sine.Amount.Alpha(c)./(Sine.Amount.Alpha(c)+Sine.Amount.NoAlpha(c));
%     cycLabel{1,c} = [num2str(round(Sine.cfg.simParams.cycles(c),2)), [' ' num2str(round(Abundance(c),2)), '']];
% end
% Figure3B.Sine_cycLabel = cellfun(@(x) strrep(x,' ','\newline'), cycLabel,'UniformOutput',false);
% 
% SNR = squeeze(Sine.SignalDetection.Amp)./squeeze(Sine.SignalDetection.fitBG);
% empiricalSNR = round(max(SNR,[],2),0);
% for a = 1:numel(Sine.cfg.simParams.amplitude)
%     Figure3B.Sine_ampSNR{a} = [num2str(Sine.cfg.simParams.amplitude(a)), ' (', num2str(empiricalSNR(a)),') '];
% end
% 
% for c = 1:numel(ExtendedB.cfg.simParams.cycles)
%     Abundance(1,c) = ExtendedB.Amount.Alpha(c)./(ExtendedB.Amount.Alpha(c)+ExtendedB.Amount.NoAlpha(c));
%     cycLabel{1,c} = [num2str(round(ExtendedB.cfg.simParams.cycles(c),2)), [' ' num2str(round(Abundance(c),2)), '']];
% end
% Figure3B.Signal_cycLabel = cellfun(@(x) strrep(x,' ','\newline'), cycLabel,'UniformOutput',false);
% 
% SNR = squeeze(ExtendedB.SignalDetection.Amp)./squeeze(ExtendedB.SignalDetection.fitBG);
% empiricalSNR = round(max(SNR,[],2),0);
% for a = 1:numel(ExtendedB.cfg.simParams.amplitude)
%     Figure3B.Signal_ampSNR{a} = [num2str(ExtendedB.cfg.simParams.amplitude(a)), ' (', num2str(empiricalSNR(a)),') '];
% end
% 
% %% save Figure data
% 
% save('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/F3B.mat', 'Figure3B');

%% load Figure data

load('/Users/kosciessa/Desktop/eBOSC/figureData/F2B.mat', 'Figure3B');

%% plot Figure

h = figure('units','normalized','position',[0 0 .6 .7]);
subplot(2,4,1); imagesc(Figure3B.Sine_oAmp); title({'Pure Sine:'; 'Overall amplitude'}); colorbar('location', 'SouthOutside');
ylabel('Simulated amplitudes'); xlabel({'Simulated cycles';'(Simulated abundance)'});
set(gca, 'YTick', 1:3:16); set(gca, 'XTick', 1:3:15); set(gca, 'YTickLabels', Figure3B.Sine_ampSNR(2:3:end)); set(gca, 'XTickLabels', Figure3B.Sine_cycLabel(1:3:14));
subplot(2,4,2); imagesc(Figure3B.Sine_Amp); title({'Pure Sine:'; 'Rhythmic amplitude'}); colorbar('location', 'SouthOutside');
ylabel('Simulated amplitudes'); xlabel({'Simulated cycles';'(Simulated abundance)'});
set(gca, 'YTick', 1:3:16); set(gca, 'XTick', 1:3:15); set(gca, 'YTickLabels', Figure3B.Sine_ampSNR(2:3:end)); set(gca, 'XTickLabels', Figure3B.Sine_cycLabel(1:3:14));
subplot(2,4,3); imagesc(Figure3B.Sine_Abn, [0 1]); title({'Pure Sine:'; 'Rhythmic abundance'}); colorbar('location', 'SouthOutside');
ylabel('Simulated amplitudes'); xlabel({'Simulated cycles';'(Simulated abundance)'});
set(gca, 'YTick', 1:3:16); set(gca, 'XTick', 1:3:15); set(gca, 'YTickLabels', Figure3B.Sine_ampSNR(2:3:end)); set(gca, 'XTickLabels', Figure3B.Sine_cycLabel(1:3:14));
subplot(2,4,4); imagesc(Figure3B.Sine_Ratio, [0 1]); title({'Pure Sine:'; 'Ratio: overall/rhythmic'}); colorbar('location', 'SouthOutside');
ylabel('Simulated amplitudes'); xlabel({'Simulated cycles';'(Simulated abundance)'});
set(gca, 'YTick', 1:3:16); set(gca, 'XTick', 1:3:15); set(gca, 'YTickLabels', Figure3B.Sine_ampSNR(2:3:end)); set(gca, 'XTickLabels', Figure3B.Sine_cycLabel(1:3:14));

% plot the same data with simulated background

subplot(2,4,5); imagesc(Figure3B.Signal_oAmp); title({'Sine+1/f:'; 'Overall amplitude'}); colorbar('location', 'SouthOutside');
ylabel('Simulated amplitudes'); xlabel({'Simulated cycles';'(Simulated abundance)'});
set(gca, 'YTick', 1:2:8); set(gca, 'XTick', 1:2:8); set(gca, 'YTickLabels', Figure3B.Signal_ampSNR(2:2:end)); set(gca, 'XTickLabels', Figure3B.Signal_cycLabel(2:2:end));
subplot(2,4,6); imagesc(Figure3B.Signal_Amp); title({'Sine+1/f:'; 'Rhythmic amplitude'}); colorbar('location', 'SouthOutside');
ylabel('Simulated amplitudes'); xlabel({'Simulated cycles';'(Simulated abundance)'});
set(gca, 'YTick', 1:2:8); set(gca, 'XTick', 1:2:8); set(gca, 'YTickLabels', Figure3B.Signal_ampSNR(2:2:end)); set(gca, 'XTickLabels', Figure3B.Signal_cycLabel(2:2:end));
subplot(2,4,7); imagesc(Figure3B.Signal_Abn, [0 1]); title({'Sine+1/f:'; 'Rhythmic abundance'}); colorbar('location', 'SouthOutside');
ylabel('Simulated amplitudes'); xlabel({'Simulated cycles';'(Simulated abundance)'});
set(gca, 'YTick', 1:2:8); set(gca, 'XTick', 1:2:8); set(gca, 'YTickLabels', Figure3B.Signal_ampSNR(2:2:end)); set(gca, 'XTickLabels', Figure3B.Signal_cycLabel(2:2:end));
subplot(2,4,8); imagesc(Figure3B.Signal_Ratio, [0 1]); title({'Sine+1/f:'; 'Ratio: overall/rhythmic'}); colorbar('location', 'SouthOutside');
ylabel('Simulated amplitudes'); xlabel({'Simulated cycles';'(Simulated abundance)'});
set(gca, 'YTick', 1:2:8); set(gca, 'XTick', 1:2:8); set(gca, 'YTickLabels', Figure3B.Signal_ampSNR(2:2:end)); set(gca, 'XTickLabels', Figure3B.Signal_cycLabel(2:2:end));

set(findall(gcf,'-property','FontSize'),'FontSize',19)

% add colorbrewer
addpath('/Volumes/LNDG/Projects/StateSwitch/dynamic/data/eeg/rest/B_analyses/A_MSE_CSD_multiVariant/T_tools/brewermap')
cBrew = brewermap(500,'RdBu');
cBrew = flipud(cBrew);
colormap(cBrew)

set(gcf,'renderer','opengl');

pn.plotFolder = '/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F_Figures/';
figureName = 'F3B';

saveas(h, [pn.plotFolder, figureName], 'fig');
saveas(h, [pn.plotFolder, figureName], 'epsc');
saveas(h, [pn.plotFolder, figureName], 'png');
