%% Plot additional simulation results for Supplemental Materials

% original: /Volumes/EEG/BOSC_Simulation/17_JQK_RED/A_scripts/G_SupplementFigure_v2.m
% similar to /Volumes/EEG/BOSC_Simulation/17_JQK_RED/A_scripts/G_ResultsPlots_extVSstnd_Paper_v10_171023.m

% %% load measures etc.
% 
% pn.data = '/Volumes/EEG/BOSC_Simulation/17_JQK_RED/D_data/D_simulation/';
% 
% Standard = load([pn.data, 'REDSimulation_standardBOSC_170630.mat']);
% ExtendedA = load([pn.data, 'REDSimulation_171023_v10A.mat']);
% ExtendedB = load([pn.data, 'REDSimulation_171023_v10B.mat']);
% ExtendedC = load([pn.data, 'REDSimulation_171023_v10C.mat']);
% ExtendedD = load([pn.data, 'REDSimulation_171023_v10D.mat']);
% 
% abundance_ep(1,:,:,:) = ExtendedA.abundance_ep(:,:,:); % rhythmic episode abundance
% abundance_ep(2,:,:,:) = ExtendedB.abundance_ep(:,:,:);
% abundance_ep(3,:,:,:) = ExtendedC.abundance_ep(:,:,:);
% abundance_ep(4,:,:,:) = ExtendedD.abundance_ep(:,:,:);
% 
% abundance_PepMeanAlpha = Standard.abundance_PepMeanAlpha; % pepisode of average detected in alpha range
% abundance_meanPep = Standard.abundance_meanPep; % average of pepisodes in alpha range
% 
% % Load information about simulation
% AmountExt = ExtendedA.Amount;
% AmountStd = Standard.Amount;
% Amount = AmountExt;
% amplitude = [0 2 4 6 8 12 16 24];
% cycles = [2 4 8 16 32 64 128 200];
% alphaFreq = 10;
% amountAlpha = round(round((cycles/alphaFreq),3)/0.004,0);
% amountAlpha(end) = 3500; % final abundance is 1, i.e. covering the entire period
% 
% % method x amplitude x cycles x repetitions [for extended BOSC]
% 
% % correct for wrong indexing on 170302:
% % correct size is 1x8x8x200
% 
% Standard.SignalDetection.Hits = Standard.SignalDetection.Hits(1,1:8,1:8,:);
% Standard.SignalDetection.Misses = Standard.SignalDetection.Misses(1,1:8,1:8,:);
% Standard.SignalDetection.CRs = Standard.SignalDetection.CRs(1,1:8,1:8,:);
% Standard.SignalDetection.FAs = Standard.SignalDetection.FAs(1,1:8,1:8,:);
% 
% abundance_PepMeanAlpha = reshape(abundance_PepMeanAlpha, [1, size(abundance_PepMeanAlpha)]);
% abundance_meanPep = reshape(abundance_meanPep, [1, size(abundance_meanPep)]);
% 
% % create abundance and cycle labels
% 
% for c = 1:numel(cycles)
%     tmp_amountAlpha = Amount.Alpha(c);
%     tmp_amountNoAlpha = Amount.NoAlpha(c);
%     Abundance(1,c) = Amount.Alpha(c)./3500;
%     cycLabel{1,c} = [num2str(round(cycles(c),2)), [' ' num2str(round(Abundance(c),2)), '']];
% end
% FigureS1.cycLabel = cellfun(@(x) strrep(x,' ','\newline'), cycLabel,'UniformOutput',false);
% 
% % create labels containing the approximated empirical SNR (overall and episode)
% % Note that the SNR will vary depending on the abundance, but this is not reflected here.
% 
% SNR = squeeze(Standard.SignalDetection.Amp)./squeeze(Standard.SignalDetection.fitBG);
% empiricalSNR = round(max(SNR,[],2),0);
% for a = 1:numel(amplitude)
%     FigureS1.ampSNR{a} = [num2str(amplitude(a)), ' (', num2str(empiricalSNR(a)),') '];
% end
% 
% FigureS1.A2 = squeeze(nanmean(Standard.SignalDetection.HitRate(1,:,:,:),4));
% FigureS1.A5 = squeeze(nanmean(ExtendedA.SignalDetection.HitRate(:,:,:),3));
% FigureS1.A8 = squeeze(nanmean(ExtendedB.SignalDetection.HitRate(:,:,:),3));
% FigureS1.A3 = squeeze(nanmean(Standard.SignalDetection.FARate(1,:,:,:),4));
% FigureS1.A6 = squeeze(nanmean(ExtendedA.SignalDetection.FARate(:,:,:),3));
% FigureS1.A9 = squeeze(nanmean(ExtendedB.SignalDetection.FARate(:,:,:),3));
% FigureS1.A1 = squeeze(nanmean(abundance_PepMeanAlpha(1,:,:,:),4))-repmat((Amount.Alpha./3500),8,1);
% FigureS1.A4 = squeeze(nanmean(abundance_ep(1,:,:,:),4))-repmat((Amount.Alpha./3500),8,1);
% FigureS1.A7 = squeeze(nanmean(abundance_ep(2,:,:,:),4))-repmat((Amount.Alpha./3500),8,1);
% FigureS1.A11 = squeeze(nanmean(ExtendedC.SignalDetection.HitRate(:,:,:),3));
% FigureS1.A14 = squeeze(nanmean(ExtendedD.SignalDetection.HitRate(:,:,:),3));
% FigureS1.A12 = squeeze(nanmean(ExtendedC.SignalDetection.FARate(:,:,:),3));
% FigureS1.A15 = squeeze(nanmean(ExtendedD.SignalDetection.FARate(:,:,:),3));
% FigureS1.A10 = squeeze(nanmean(abundance_ep(3,:,:,:),4))-repmat((Amount.Alpha./3500),8,1);
% FigureS1.A13 = squeeze(nanmean(abundance_ep(4,:,:,:),4))-repmat((Amount.Alpha./3500),8,1);
% 
% %% save Figure data
% 
% save('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/S1.mat', 'FigureS1');

%% load Figure data

load('/Users/kosciessa/Desktop/eBOSC/figureData/S2.mat', 'FigureS1');

addpath(genpath('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F4_FigureTools/'))

%% create overview plot

h = figure('units','normalized','position',[0 0 1 1]);
subplot(5,3,2);
    imagesc(FigureS1.A2, [0 1]);
    addNumericalValues_Figure(FigureS1.A2, FigureS1.ampSNR, FigureS1.cycLabel)
    title('Hit Rate: Standard BOSC');
subplot(5,3,5);
    imagesc(FigureS1.A5, [0 1]);
    addNumericalValues_Figure(FigureS1.A5, FigureS1.ampSNR, FigureS1.cycLabel)
    title('Hit Rate: MaxBias raw');
subplot(5,3,8);
    imagesc(FigureS1.A8, [0 1]);
    addNumericalValues_Figure(FigureS1.A8, FigureS1.ampSNR, FigureS1.cycLabel)
    title('Hit Rate: MaxBias PT');
subplot(5,3,3);
    imagesc(FigureS1.A3, [0 1]);
    addNumericalValues_Figure(FigureS1.A3, FigureS1.ampSNR, FigureS1.cycLabel)
    title('FA Rate: Standard BOSC');
subplot(5,3,6);
    imagesc(FigureS1.A6, [-.2 .2]);
    addNumericalValues_Figure(FigureS1.A6, FigureS1.ampSNR, FigureS1.cycLabel)
    title('FA Rate: MaxBias raw');
subplot(5,3,9);
    imagesc(FigureS1.A9, [-.2 .2]);
    addNumericalValues_Figure(FigureS1.A9, FigureS1.ampSNR, FigureS1.cycLabel)
    title('FA Rate: MaxBias PT');
subplot(5,3,1);
    imagesc(FigureS1.A1, [-.2, .2]);
    addNumericalValues_Figure(FigureS1.A1, FigureS1.ampSNR, FigureS1.cycLabel)
    title('Abundance Error: Standard BOSC');
subplot(5,3,4);
    imagesc(FigureS1.A4, [-.2, .2]);
    addNumericalValues_Figure(FigureS1.A4, FigureS1.ampSNR, FigureS1.cycLabel)
    title('Abundance Error: MaxBias raw');
subplot(5,3,7);
    imagesc(FigureS1.A7, [-.2, .2]);
    addNumericalValues_Figure(FigureS1.A7, FigureS1.ampSNR, FigureS1.cycLabel)
    title('Abundance Error: MaxBias PT');
subplot(5,3,11);
    imagesc(FigureS1.A11, [0 1]);
    addNumericalValues_Figure(FigureS1.A11, FigureS1.ampSNR, FigureS1.cycLabel)
    title('Hit Rate: FWHM raw');
subplot(5,3,14);
    imagesc(FigureS1.A14, [0 1]);
    addNumericalValues_Figure(FigureS1.A14, FigureS1.ampSNR, FigureS1.cycLabel)
    title('Hit Rate: FWHM PT');
subplot(5,3,12);
    imagesc(FigureS1.A12, [-.2 .2]);
    addNumericalValues_Figure(FigureS1.A12, FigureS1.ampSNR, FigureS1.cycLabel)
    title('FA Rate: FWHM raw');
subplot(5,3,15);
    imagesc(FigureS1.A15, [-.2 .2]);
    addNumericalValues_Figure(FigureS1.A15, FigureS1.ampSNR, FigureS1.cycLabel)
    title('FA Rate: FWHM PT');
subplot(5,3,10);
    imagesc(FigureS1.A10, [-.2, .2]);
    addNumericalValues_Figure(FigureS1.A10, FigureS1.ampSNR, FigureS1.cycLabel)
    title('Abundance Error: FWHM raw');
subplot(5,3,13);
    imagesc(FigureS1.A13, [-.2, .2]);
    addNumericalValues_Figure(FigureS1.A13, FigureS1.ampSNR, FigureS1.cycLabel)
    title('Abundance Error: FWHM PT');

set(findall(gcf,'-property','FontSize'),'FontSize',12)

%% save Figure

pn.plotFolder = '/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F_Figures/';
figureName = 'S1';

saveas(h, [pn.plotFolder, figureName], 'fig');
saveas(h, [pn.plotFolder, figureName], 'epsc');
saveas(h, [pn.plotFolder, figureName], 'png');