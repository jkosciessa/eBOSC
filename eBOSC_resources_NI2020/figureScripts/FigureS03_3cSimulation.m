%% plot results from simulation with 3 cycles
% previously /Volumes/EEG/BOSC_Simulation/17_JQK_RED/A_scripts/K_simulation3c/A_scripts/B_Plot_detectionPerformance_3c.m

% %% load measures etc.
% 
% pn.data = '/Volumes/EEG/BOSC_Simulation/17_JQK_RED/A_scripts/K_simulation3c/B_data/';
% pn.plotFolder = '/Volumes/EEG/BOSC_Simulation/17_JQK_RED/A_scripts/K_simulation3c/C_figures/';
% 
% % Load simulation results
% 
% ExtendedB = load([pn.data, 'REDSimulation_181005_B.mat']);
% 
% % create abundance and cycle labels & labels containing the approximated empirical SNR (overall and episode)
% % Note that the SNR will vary depending on the abundance, but this is not reflected here.
% 
% Amount = ExtendedB.Amount;
% FigureS6.amplitude = [0 2 4 6 8 12 16 24];
% cycles = [2 4 8 16 32 64 128 200];
% 
% for c = 1:numel(cycles)
%     tmp_amountAlpha = Amount.Alpha(c);
%     tmp_amountNoAlpha = Amount.NoAlpha(c);
%     Abundance(1,c) = Amount.Alpha(c)./3500;
%     cycLabel{1,c} = [num2str(round(cycles(c),2)), [' ' num2str(round(Abundance(c),2)), '']];
% end
% FigureS6.cycLabel = cellfun(@(x) strrep(x,' ','\newline'), cycLabel,'UniformOutput',false);
% 
% SNR = squeeze(ExtendedB.SignalDetection.Amp)./squeeze(ExtendedB.SignalDetection.fitBG);
% empiricalSNR = round(max(SNR,[],2),0);
% for a = 1:numel(amplitude)
%     FigureS6.ampSNR{a} = [num2str(amplitude(a)), ' (', num2str(empiricalSNR(a)),') '];
% end
% 
% %% prepare Figure data
% 
% FigureS6.A1 = squeeze(nanmean(ExtendedB.abundance_ep(:,:,:),3))-repmat((Amount.Alpha./3500),8,1);
% FigureS6.A2 = squeeze(nanmean(ExtendedB.SignalDetection.HitRate(:,:,:),3));
% FigureS6.A3 = squeeze(nanmean(ExtendedB.SignalDetection.FARate(:,:,:),3));
% 
% %% save Figure data
% 
% save('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/S6.mat', 'FigureS6')

%% load Figure data

load('/Users/kosciessa/Desktop/eBOSC/figureData/S3.mat', 'FigureS6')

%% Figure SXX

h = figure('units','normalized','position',[0 0 1 .6]);
subplot(2,3,1);
    imagesc(FigureS6.A1, [-.2, .2]);
    addNumericalValues_Figure(FigureS6.A1, FigureS6.amplitude, FigureS6.cycLabel)
    title('Abundance Error: Extended BOSC 3 cycle Wavelet');
subplot(2,3,2);
    imagesc(FigureS6.A2, [0 1]);
    addNumericalValues_Figure(FigureS6.A2, FigureS6.amplitude, FigureS6.cycLabel)
    title('Hit Rate: Extended BOSC 3 cycle Wavelet');
subplot(2,3,3);
    imagesc(FigureS6.A3, [-.2 .2]);
    addNumericalValues_Figure(FigureS6.A3, FigureS6.amplitude, FigureS6.cycLabel)
    title('FA Rate: Extended BOSC 3 cycle Wavelet');
    
%% load 7 cycle results and plot comparison

load('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/F4.mat', 'Figure4')

subplot(2,3,4);
    imagesc(FigureS6.A1-Figure4.A12, [-1 1]);
    addNumericalValues_Figure(FigureS6.A1-Figure4.A12, FigureS6.amplitude, FigureS6.cycLabel)
    title('Abundance Error: Extended BOSC 3c vs. 6c');
subplot(2,3,5);
    imagesc(FigureS6.A2-Figure4.A22, [-1 1]);
    addNumericalValues_Figure(FigureS6.A2-Figure4.A22, FigureS6.amplitude, FigureS6.cycLabel)
    title('Hit Rate: Extended BOSC 3c vs. 6c');
subplot(2,3,6);
    imagesc(FigureS6.A3-Figure4.A32, [-1 1]);
    addNumericalValues_Figure(FigureS6.A3-Figure4.A32, FigureS6.amplitude, FigureS6.cycLabel)
    title('FA Rate: Extended BOSC 3c vs. 6c');

set(findall(gcf,'-property','FontSize'),'FontSize',15)

pn.plotFolder = '/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F_Figures/';
figureName = 'S6';

saveas(h, [pn.plotFolder, figureName], 'fig');
saveas(h, [pn.plotFolder, figureName], 'epsc');
saveas(h, [pn.plotFolder, figureName], 'png');
