%% Figure 4: plot simulation performance of BOSC and eBOSC

% %% load measures etc.
% 
% pn.data = '/Volumes/EEG/BOSC_Simulation/17_JQK_RED/D_data/D_simulation/';
% pn.plotFolder = '/Volumes/EEG/BOSC_Simulation/17_JQK_RED/B_figures/';
% 
% % Load simulation results
% 
% Standard = load([pn.data, 'REDSimulation_standardBOSC_170630.mat']);
% ExtendedB = load([pn.data, 'REDSimulation_171023_v10B.mat']);
% 
% % correct for wrong indexing: correct size is 1x8x8x200
% 
% Standard.SignalDetection.Hits = Standard.SignalDetection.Hits(1,1:8,1:8,:);
% Standard.SignalDetection.Misses = Standard.SignalDetection.Misses(1,1:8,1:8,:);
% Standard.SignalDetection.CRs = Standard.SignalDetection.CRs(1,1:8,1:8,:);
% Standard.SignalDetection.FAs = Standard.SignalDetection.FAs(1,1:8,1:8,:);
% 
% % create abundance and cycle labels & labels containing the approximated empirical SNR (overall and episode)
% % Note that the SNR will vary depending on the abundance, but this is not reflected here.
% 
% Amount = ExtendedB.Amount;
% amplitude = [0 2 4 6 8 12 16 24];
% cycles = [2 4 8 16 32 64 128 200];
% 
% for c = 1:numel(cycles)
%     tmp_amountAlpha = Amount.Alpha(c);
%     tmp_amountNoAlpha = Amount.NoAlpha(c);
%     Abundance(1,c) = Amount.Alpha(c)./3500;
%     cycLabel{1,c} = [num2str(round(cycles(c),2)), [' ' num2str(round(Abundance(c),2)), '']];
% end
% Figure4.cycLabel = cellfun(@(x) strrep(x,' ','\newline'), cycLabel,'UniformOutput',false);
% 
% SNR = squeeze(Standard.SignalDetection.Amp)./squeeze(Standard.SignalDetection.fitBG);
% empiricalSNR = round(max(SNR,[],2),0);
% for a = 1:numel(amplitude)
%     Figure4.ampSNR{a} = [num2str(amplitude(a)), ' (', num2str(empiricalSNR(a)),') '];
% end
% 
% %% prepare Figure data
% 
% Figure4.A11 = squeeze(nanmean(Standard.abundance_PepMeanAlpha(:,:,:),3))-repmat((Amount.Alpha./3500),8,1);
% Figure4.A12 = squeeze(nanmean(ExtendedB.abundance_ep(:,:,:),3))-repmat((Amount.Alpha./3500),8,1);
% Figure4.A21 = squeeze(nanmean(Standard.SignalDetection.HitRate(1,:,:,:),4));
% Figure4.A22 = squeeze(nanmean(ExtendedB.SignalDetection.HitRate(:,:,:),3));
% Figure4.A31 = squeeze(nanmean(Standard.SignalDetection.FARate(1,:,:,:),4));
% Figure4.A32 = squeeze(nanmean(ExtendedB.SignalDetection.FARate(:,:,:),3));
% 
% Figure4.C.Standard_fitBG = Standard.SignalDetection.fitBG(:);
% Figure4.C.Standard_oAmp = Standard.SignalDetection.oAmp(:);
% Figure4.C.ExtendedB_fitBG = ExtendedB.SignalDetection.fitBG(:);
% Figure4.C.ExtendedB_oAmp = ExtendedB.SignalDetection.oAmp(:);
% 
% Figure4.D.Standard_abn = Standard.SignalDetection.abn(:);
% Figure4.D.Standard_oAmp = Standard.SignalDetection.oAmp(:);
% Figure4.D.ExtendedB_abn = ExtendedB.SignalDetection.abn(:);
% Figure4.D.ExtendedB_oAmp = ExtendedB.SignalDetection.oAmp(:);
% 
% % calculate between-trial amplitude abundance associations
% for a = 1:8
%     for c = 1:8
%         % for Standard BOSC
%         signalOfInterest_Amp = squeeze(Standard.SignalDetection.oAmp(1,a,c,:));
%         signalOfInterest_Abn = squeeze(Standard.SignalDetection.abn(1,a,c,:));
%         r = corrcoef(signalOfInterest_Amp(~isnan(signalOfInterest_Amp)), signalOfInterest_Abn(~isnan(signalOfInterest_Amp)));
%         if isnan(r(2))
%             r(2) = 0; % if one variable has no variance, set the correlation to zero
%         end
%         Figure4.E1.ZMat(1,a,c) = .5*log((1+r(2))./(1-r(2))); % Fisher's z-transform
%         % for eBOSC
%         signalOfInterest_Amp = squeeze(ExtendedB.SignalDetection.oAmp(a,c,:));
%         signalOfInterest_Abn = squeeze(ExtendedB.SignalDetection.abn(a,c,:));
%         r = corrcoef(signalOfInterest_Amp(~isnan(signalOfInterest_Amp)), signalOfInterest_Abn(~isnan(signalOfInterest_Amp)));
%         if isnan(r(2))
%             r(2) = 0; % if one variable has no variance, set the correlation to zero
%         end
%         Figure4.E1.ZMat(2,a,c) = .5*log((1+r(2))./(1-r(2))); % Fisher's z-transform
%     end
% end
% 
% % condition: 1-pure sine; 2- sine with BG
% 
% Figure4.B.colorm = brewermap(4,'RdYlBu'); % use external function to create colormap
% 
% AmpAbnSimulation = load(['/Volumes/EEG/BOSC_Simulation/17_JQK_RED/A_scripts/I_amplitudeRecovery/B_data/AmpAbnSimulation_v5.mat']);
% 
% SNR = squeeze(AmpAbnSimulation.SignalDetection.Amp(:,:,4,:))./squeeze(AmpAbnSimulation.SignalDetection.fitBG(:,:,4,:));
% empiricalSNR = round(nanmean(SNR,2),1);
% for a = 1:numel(AmpAbnSimulation.cfg.simParams.amplitude)
%     Figure4.B.ampSNR{a} = [num2str(AmpAbnSimulation.cfg.simParams.amplitude(a)), ' (', num2str(empiricalSNR(a)),') '];
% end
% 
% SNRo = squeeze(AmpAbnSimulation.SignalDetection.oAmp(:,:,4,:))./squeeze(AmpAbnSimulation.SignalDetection.oBGestimate(:,:,4,:));
% empiricalSNRo = round(nanmean(SNRo(:,:),2),1);
% for a = 1:numel(AmpAbnSimulation.cfg.simParams.amplitude)
%     Figure4.B.ampSNR{a} = [num2str(AmpAbnSimulation.cfg.simParams.amplitude(a)), ' (', num2str(empiricalSNRo(a)),') '];
% end
% 
% Figure4.B1.Sim1_Abn = mean(AmpAbnSimulation.SignalDetection.abn(:,:,4,:),4);
% Figure4.B1.Sim1_Abn_noBG = mean(AmpAbnSimulation.SignalDetection.abn(:,:,2,:),4);
% Figure4.B1.Sim1_simAmp = mean(AmpAbnSimulation.SignalDetection.Amp(:,:,2,:),4);
% Figure4.B1.Sim1_eAmp = mean(AmpAbnSimulation.SignalDetection.Amp(:,:,4,:),4);
% Figure4.B1.Sim1_oAmp = mean(AmpAbnSimulation.SignalDetection.oAmp(:,:,4,:),4);
% Figure4.B1.Sim1_fitBG = mean(AmpAbnSimulation.SignalDetection.fitBG(:,:,4,:),4);
% 
% Figure4.B2.Sim25_Abn = mean(AmpAbnSimulation.SignalDetection.abn(:,:,3,:),4);
% Figure4.B2.Sim25_Abn_noBG = mean(AmpAbnSimulation.SignalDetection.abn(:,:,1,:),4);
% Figure4.B2.Sim25_simAmp = mean(AmpAbnSimulation.SignalDetection.Amp(:,:,1,:),4);
% Figure4.B2.Sim25_eAmp = mean(AmpAbnSimulation.SignalDetection.Amp(:,:,3,:),4);
% Figure4.B2.Sim25_oAmp = mean(AmpAbnSimulation.SignalDetection.oAmp(:,:,3,:),4);
% Figure4.B2.Sim25_fitBG = mean(AmpAbnSimulation.SignalDetection.fitBG(:,:,3,:),4);
% 
% Figure4.B2.SimAbn = mean(AmpAbnSimulation.SignalDetection.abn(:,:,1,:),4);
% Figure4.B2.SimAmp = mean(AmpAbnSimulation.SignalDetection.Amp(:,:,1,:),4);
% 
% % save amp-abn relationship for comparison with empirical resting state
% 
% % % figure;
% % % plot(mean(AmpAbnSimulation.SignalDetection.abn(:,:,4,:),4), empiricalSNR);
% % 
% % simData = [];
% % simData.Abn = mean(SignalDetection.abn(:,:,4,:),4);
% % simData.SNR = round(max(SNR,[],2),0);
% % 
% % save([pn.data, '_SimDataforRestingStateComparison_v4.mat'], 'simData');
% % 
% % % figure;
% % % plot(mean(AmpAbnSimulation.SignalDetection.abn(:,:,4,:),4), empiricalSNRo);
% % 
% % simData = [];
% % simData.Abn = mean(SignalDetection.abn(:,:,4,:),4);
% % simData.SNR = empiricalSNRo;
% % 
% % save([pn.data, '_SimDataforRestingStateComparison_overall_v4.mat'], 'simData');
% 
% %% save Figure data
% 
% save('/Volumes/EEG/BOSC/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/F4.mat', 'Figure4')

%% load Figure data

load('/Users/kosciessa/Desktop/eBOSC/figureData/F3.mat', 'Figure4')

addpath(genpath('/Volumes/EEG/BOSC/BOSC_SternRest/X_documentation/B_2018_Manuscript/F4_FigureTools/'))

%% Figure A

h = figure('units','normalized','position',[0 0 1 1]);
subplot(3,4,1);
    imagesc(Figure4.A11, [-.2, .2]);
    %ax = gca; colorData = ax.Children.CData;
    addNumericalValues_Figure(Figure4.A11, Figure4.ampSNR, Figure4.cycLabel)
    title('Abundance Error: Standard BOSC');
subplot(3,4,2);
    imagesc(Figure4.A12, [-.2, .2]);
    addNumericalValues_Figure(Figure4.A12, Figure4.ampSNR, Figure4.cycLabel)
    title('Abundance Error: Extended BOSC');
subplot(3,4,4+1);
    imagesc(Figure4.A21, [0 1]);
    addNumericalValues_Figure(Figure4.A21, Figure4.ampSNR, Figure4.cycLabel)
    title('Hit Rate: Standard BOSC');
subplot(3,4,4+2);
    imagesc(Figure4.A22, [0 1]);
    addNumericalValues_Figure(Figure4.A22, Figure4.ampSNR, Figure4.cycLabel)
    title('Hit Rate: Extended BOSC');
subplot(3,4,8+1);
    imagesc(Figure4.A31, [0 1]);
    addNumericalValues_Figure(Figure4.A31, Figure4.ampSNR, Figure4.cycLabel)
    title('FA Rate: Standard BOSC');
subplot(3,4,8+2);
    imagesc(Figure4.A32, [-.2 .2]);
    addNumericalValues_Figure(Figure4.A32, Figure4.ampSNR, Figure4.cycLabel)
    title('FA Rate: Extended BOSC');

%% Figure 3 C, D: oAmp-fitBG relationship across all simulated amplitudes and cycles

subplot(3,4,4+3);
    hold on;
    h1 = scatter(Figure4.C.Standard_fitBG, Figure4.C.Standard_oAmp, 8, 'k', 'filled'); xlabel('Estimated Background amplitude'); ylabel('Overall Alpha amplitude'); xlim([21 35])
    h2 = scatter(Figure4.C.ExtendedB_fitBG, Figure4.C.ExtendedB_oAmp, 8, 'r', 'filled'); xlabel('Estimated Background amplitude'); ylabel('Overall Alpha amplitude'); xlim([21 35])
    [r1, p1] = corrcoef(Figure4.C.Standard_fitBG(~isnan(Figure4.C.Standard_fitBG)), Figure4.C.Standard_oAmp(~isnan(Figure4.C.Standard_fitBG)));
    [r2, p2] = corrcoef(Figure4.C.ExtendedB_fitBG(~isnan(Figure4.C.ExtendedB_fitBG)), Figure4.C.ExtendedB_oAmp(~isnan(Figure4.C.ExtendedB_fitBG)));
    legend([h1, h2], {['Standard BOSC: r = ',num2str(round(r1(2),2)),'; p <.001']; ['extended BOSC: r = ',num2str(round(r2(2),2)),'; p <.001']}, 'location', 'NorthWest'); legend('boxoff')
    ylim([0 275]);
    title({'eBOSC decreases rhythmic bias';'on background estimates'})
subplot(3,4,4+4);
    hold on;
    h1 = scatter(Figure4.D.Standard_oAmp(:), Figure4.D.Standard_abn(:), 8, 'k', 'filled'); xlabel('Overall Alpha amplitude'); ylabel('Estimated Abundance'); % xlim([21 35])
    h2 = scatter(Figure4.D.ExtendedB_oAmp(:), Figure4.D.ExtendedB_abn(:), 8, 'r', 'filled'); xlabel('Overall Alpha amplitude'); ylabel('Estimated Abundance');% xlim([21 35])
    line([20 220],[.01, .01], 'Color', 'k')
    line([20 220],[.03, .03], 'Color', 'k')
    line([20 220],[.06, .06], 'Color', 'k')
    line([20 220],[.11, .11], 'Color', 'k')
    line([20 220],[.23, .23], 'Color', 'k')
    line([20 220],[.46, .46], 'Color', 'k')
    line([20 220],[.91, .91], 'Color', 'k')
    line([20 220],[1, 1], 'Color', 'k')
    xlim([20 220])
    title({'eBOSC more accurately captures';'simulated durations'})
 
%% Figure 3 E: correlations within simulated amplitude and abundance

subplot(3,4,8+3);
    imagesc(squeeze(Figure4.E1.ZMat(1,:,:)),[-1 1]);
    addNumericalValues_Figure(squeeze(Figure4.E1.ZMat(1,:,:)), Figure4.ampSNR, Figure4.cycLabel)
    title({'Standard BOSC:';'Trial-wise amplitude-abundance association'});
subplot(3,4,8+4);
    imagesc(squeeze(Figure4.E1.ZMat(2,:,:)),[-1 1]);
    addNumericalValues_Figure(squeeze(Figure4.E1.ZMat(2,:,:)), Figure4.ampSNR, Figure4.cycLabel)
    title({'Extended BOSC:';'Trial-wise amplitude-abundance association'});

%% Figure 1B: simulated amplitude and abundance

subplot(3,4,3);
    l1 = scatter(Figure4.B1.Sim1_Abn, Figure4.B1.Sim1_simAmp,60,'MarkerEdgeColor', [1,1,1], 'MarkerFaceColor', [0,0,0]); hold on; % simulated amplitudes
    l2 = scatter(Figure4.B1.Sim1_Abn_noBG, Figure4.B1.Sim1_simAmp,60,'MarkerEdgeColor', [1,1,1], 'MarkerFaceColor', [.5 .5 .5 ]); hold on; % simulated amplitudes
    l3 = scatter(Figure4.B1.Sim1_Abn, Figure4.B1.Sim1_oAmp,60,'MarkerEdgeColor', [1,1,1], 'MarkerFaceColor', Figure4.B.colorm(4,:)); hold on; % overall amplitudes
    l4 = scatter(Figure4.B1.Sim1_Abn, Figure4.B1.Sim1_eAmp,60, 'MarkerEdgeColor', [1,1,1], 'MarkerFaceColor', Figure4.B.colorm(1,:)); % detected rhythmic amplitudes
    l5 = scatter(Figure4.B1.Sim1_Abn, Figure4.B1.Sim1_eAmp-Figure4.B1.Sim1_fitBG,60, 'MarkerEdgeColor', [1,1,1], 'MarkerFaceColor', Figure4.B.colorm(2,:)); % detected rhythmic amplitude-BG
    ylabel('Empirical Amplitude [a.u.]'); xlabel('Empirical Abundance');
    title({'Simulated vs. empirical amplitudes';'Simulated Abundance: 1'});
    line([1 1], [0, 150], 'Color','k', 'LineWidth', 2) % add simulated abundance
    ylim([25 150])
    legend([l2, l1, l3], {'simulated data'; ['simulated data',char(10),'(interpolated abundance)']; 'overall estimate'}, 'Location', 'northwest')
    legend('boxoff')
    
% fit to lower simulated abundance

subplot(3,4,4); cla;
    l1 = scatter(Figure4.B2.Sim25_Abn, Figure4.B2.Sim25_simAmp,60,'MarkerEdgeColor', [1,1,1], 'MarkerFaceColor', [0,0,0]); hold on; % simulated amplitudes
    l2 = scatter(Figure4.B2.Sim25_Abn_noBG, Figure4.B2.Sim25_simAmp,60,'MarkerEdgeColor', [1,1,1], 'MarkerFaceColor', [.5 .5 .5 ]); hold on; % simulated amplitudes
    l3 = scatter(Figure4.B2.Sim25_Abn, Figure4.B2.Sim25_oAmp,60,'MarkerEdgeColor', [1,1,1], 'MarkerFaceColor', Figure4.B.colorm(4,:)); hold on; % overall amplitudes
    l4 = scatter(Figure4.B2.Sim25_Abn, Figure4.B2.Sim25_eAmp,60, 'MarkerEdgeColor', [1,1,1], 'MarkerFaceColor', Figure4.B.colorm(1,:)); % detected rhythmic amplitudes
    l5 = scatter(Figure4.B2.Sim25_Abn, Figure4.B2.Sim25_eAmp-Figure4.B2.Sim25_fitBG,60, 'MarkerEdgeColor', [1,1,1], 'MarkerFaceColor', Figure4.B.colorm(2,:)); % detected rhythmic amplitude-BG
    ylabel('Empirical Amplitude [a.u.]'); xlabel('Empirical Abundance');
    title({'Simulated vs. empirical amplitudes';'Simulated Abundance: .25'});
    line([.25 .25], [0, 150], 'Color','k', 'LineWidth', 2) % add simulated abundance 
    ylim([25 150])
    legend([l4, l5], {'rhythmic estimate'; ['rhythmic estimate' char(10) 'excl. background']}, 'Location', 'northwest')
    legend('boxoff')
    
set(findall(gcf,'-property','FontSize'),'FontSize',20)

subplot(3,4,1); set(gca, 'XTick', 1:8); set(findall(gca,'-property','FontSize'),'FontSize',17); set(gca, 'FontSize', 20);
subplot(3,4,2); set(gca, 'XTick', 1:8); set(findall(gca,'-property','FontSize'),'FontSize',17); set(gca, 'FontSize', 20);
subplot(3,4,5); set(gca, 'XTick', 1:8); set(findall(gca,'-property','FontSize'),'FontSize',17); set(gca, 'FontSize', 20);
subplot(3,4,6); set(gca, 'XTick', 1:8); set(findall(gca,'-property','FontSize'),'FontSize',17); set(gca, 'FontSize', 20);
subplot(3,4,9); set(gca, 'XTick', 1:8); set(findall(gca,'-property','FontSize'),'FontSize',17); set(gca, 'FontSize', 20);
subplot(3,4,10); set(gca, 'XTick', 1:8); set(findall(gca,'-property','FontSize'),'FontSize',17); set(gca, 'FontSize', 20);
subplot(3,4,11); set(gca, 'XTick', 1:8); set(findall(gca,'-property','FontSize'),'FontSize',17); set(gca, 'FontSize', 20);
subplot(3,4,12); set(gca, 'XTick', 1:8); set(findall(gca,'-property','FontSize'),'FontSize',17); set(gca, 'FontSize', 20);

subplot(3,4,4);
% add SNR as text
x = Figure4.B2.SimAbn;
y = Figure4.B2.SimAmp;
z = Figure4.B.ampSNR;
for K = 1:2:numel(x)
    t(K) = text(x(K)+.02,y(K),z{K}, 'FontSize', 14);
end

% change colormap
addpath('/Volumes/EEG/BOSC_Sternberg/T_tools/brewermap')
cBrew = brewermap(500,'RdBu');
cBrew = flipud(cBrew);
colormap(cBrew)

%% save Figure

pn.plotFolder = '/Volumes/EEG/BOSC/BOSC_SternRest/X_documentation/B_2018_Manuscript/Figures_resubmit1/F1_Figures/';
figureName = 'F3_v2';

saveas(h, [pn.plotFolder, figureName], 'fig');
saveas(h, [pn.plotFolder, figureName], 'epsc');
saveas(h, [pn.plotFolder, figureName], 'png');
