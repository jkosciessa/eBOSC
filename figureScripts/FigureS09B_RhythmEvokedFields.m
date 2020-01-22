%% Plot Rhythm-evoked fields at the target frequency

% load('/Users/kosciessa/Desktop/mntTardis/Stern_WIP/B_eBOSC_180509_full_wl6_noDur/B_180509_MergeRhythmCharacteristics/B_data/A5_FigureData.mat', 'FigureXX')
% 
% Figure15B.time = FigureXX.time;
% Figure15B.ThetaPowerThetaOnset = zscore(squeeze(nanmedian(nanmedian(FigureXX.ThetaOnset_diag(:,:,1,:,:),4),1)),[],2);
% Figure15B.ThetaPowerThetaOffset = zscore(squeeze(nanmedian(nanmedian(FigureXX.ThetaOffset_diag(:,:,1,:,:),4),1)),[],2);
% Figure15B.ThetaPowerThetaOnOffset = zscore(squeeze(nanmedian(nanmedian(FigureXX.ThetaOnMinOff_diag(:,:,1,:,:),4),1)),[],2);
% 
% Figure15B.AlphaPowerAlphaOnset = zscore(squeeze(nanmedian(nanmedian(FigureXX.AlphaOnset_diag(:,:,2,:,:),4),1)),[],2);
% Figure15B.AlphaPowerAlphaOffset = zscore(squeeze(nanmedian(nanmedian(FigureXX.AlphaOffset_diag(:,:,2,:,:),4),1)),[],2);
% Figure15B.AlphaPowerAlphaOnOffset = zscore(squeeze(nanmedian(nanmedian(FigureXX.AlphaOnMinOff_diag(:,:,2,:,:),4),1)),[],2);
% 
% Figure15B.BetaPowerBetaOnset = zscore(squeeze(nanmedian(nanmedian(FigureXX.BetaOnset_diag(:,:,3,:,:),4),1)),[],2);
% Figure15B.BetaPowerBetaOffset = zscore(squeeze(nanmedian(nanmedian(FigureXX.BetaOffset_diag(:,:,3,:,:),4),1)),[],2);
% Figure15B.BetaPowerBetaOnOffset = zscore(squeeze(nanmedian(nanmedian(FigureXX.BetaOnMinOff_diag(:,:,3,:,:),4),1)),[],2);
% 
% Figure15B.GammaPowerGammaOnset = zscore(squeeze(nanmedian(nanmedian(FigureXX.GammaOnset_diag(:,:,4,:,:),4),1)),[],2);
% Figure15B.GammaPowerGammaOffset = zscore(squeeze(nanmedian(nanmedian(FigureXX.GammaOffset_diag(:,:,4,:,:),4),1)),[],2);
% Figure15B.GammaPowerGammaOnOffset = zscore(squeeze(nanmedian(nanmedian(FigureXX.GammaOnMinOff_diag(:,:,4,:,:),4),1)),[],2);
% 
% %% save Figure data
% 
% save('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/F15B.mat', 'Figure15B');
% 
%% load Figure data

load('/Users/kosciessa/Desktop/eBOSC/figureData/S9B.mat', 'Figure15B');

%% plot Figure

h = figure('units','normalized','position',[0 .1 .9 .7]);

subplot(3,5,0*5+1);
    imagesc(Figure15B.time, [], Figure15B.ThetaPowerThetaOnset);
    line([0 0], [60 0], 'Color', [0 0 0]); colorbar;
    xlabel('Time centered on Onset (s)'); ylabel('Subjects'); title('Onset: Theta')
subplot(3,5,1*5+1);
    imagesc(Figure15B.time, [], Figure15B.ThetaPowerThetaOffset);
    line([0 0], [60 0], 'Color', [0 0 0]); colorbar;
    xlabel('Time centered on Offset (s)'); ylabel('Subjects'); title('Offset: Theta')
subplot(3,5,2*5+1);
    imagesc(Figure15B.time, [], Figure15B.ThetaPowerThetaOnOffset);
    line([0 0], [60 0], 'Color', [0 0 0]); colorbar;
    xlabel('Time centered on On/Offset (s)'); ylabel('Subjects'); title('On-Offset: Theta')

subplot(3,5,0*5+2);
    imagesc(Figure15B.time, [], Figure15B.AlphaPowerAlphaOnset);
    line([0 0], [60 0], 'Color', [0 0 0]); colorbar;
    xlabel('Time centered on Onset (s)'); ylabel('Subjects'); title('Onset: Alpha')
subplot(3,5,1*5+2);
    imagesc(Figure15B.time, [], Figure15B.AlphaPowerAlphaOffset);
    line([0 0], [60 0], 'Color', [0 0 0]); colorbar;
    xlabel('Time centered on Offset (s)'); ylabel('Subjects'); title('Offset: Alpha')
subplot(3,5,2*5+2);
    imagesc(Figure15B.time, [], Figure15B.AlphaPowerAlphaOnOffset);
    line([0 0], [60 0], 'Color', [0 0 0]); colorbar;
    xlabel('Time centered on On/Offset (s)'); ylabel('Subjects'); title('On-Offset: Alpha')

subplot(3,5,0*5+3);
    imagesc(Figure15B.time, [], Figure15B.BetaPowerBetaOnset);
    line([0 0], [60 0], 'Color', [0 0 0]); colorbar;
    xlabel('Time centered on Onset (s)'); ylabel('Subjects'); title('Onset: Beta')
subplot(3,5,1*5+3);
    imagesc(Figure15B.time, [], Figure15B.BetaPowerBetaOffset);
    line([0 0], [60 0], 'Color', [0 0 0]); colorbar;
    xlabel('Time centered on Offset (s)'); ylabel('Subjects'); title('Offset: Beta')
subplot(3,5,2*5+3);
    imagesc(Figure15B.time, [], Figure15B.BetaPowerBetaOnOffset);
    line([0 0], [60 0], 'Color', [0 0 0]); colorbar;
    xlabel('Time centered on On/Offset (s)'); ylabel('Subjects'); title('On-Offset: Beta')

subplot(3,5,0*5+4);
    imagesc(Figure15B.time, [], Figure15B.GammaPowerGammaOnset);
    line([0 0], [60 0], 'Color', [0 0 0]); colorbar;
    xlabel('Time centered on Onset (s)'); ylabel('Subjects'); title('Onset: Gamma')
subplot(3,5,1*5+4);
    imagesc(Figure15B.time, [], Figure15B.GammaPowerGammaOffset);
    line([0 0], [60 0], 'Color', [0 0 0]); colorbar;
    xlabel('Time centered on Offset (s)'); ylabel('Subjects'); title('Offset: Gamma')
subplot(3,5,2*5+4);
    imagesc(Figure15B.time, [], Figure15B.GammaPowerGammaOnOffset);
    line([0 0], [60 0], 'Color', [0 0 0]); colorbar;
    xlabel('Time centered on On/Offset (s)'); ylabel('Subjects'); title('On-Offset: Gamma')

subplot(3,5,5);
hold on;
plot(Figure15B.time, squeeze(nanmean(Figure15B.ThetaPowerThetaOnset,1)), 'LineWidth', 2);
plot(Figure15B.time, squeeze(nanmean(Figure15B.AlphaPowerAlphaOnset,1)), 'LineWidth', 2);
plot(Figure15B.time, squeeze(nanmean(Figure15B.BetaPowerBetaOnset,1)), 'LineWidth', 2);
plot(Figure15B.time, squeeze(nanmean(Figure15B.GammaPowerGammaOnset,1)), 'LineWidth', 2);
line([0 0], [-4.5 4.5], 'Color', [0 0 0]); xlim([-1.2 1.2]); ylim([-4.5 4.5]);
xlabel('Time centered on Onset (s)'); ylabel('Power (z-score)'); title('Onset')

subplot(3,5,10);
hold on;
plot(Figure15B.time, squeeze(nanmean(Figure15B.ThetaPowerThetaOffset,1)), 'LineWidth', 2);
plot(Figure15B.time, squeeze(nanmean(Figure15B.AlphaPowerAlphaOffset,1)), 'LineWidth', 2);
plot(Figure15B.time, squeeze(nanmean(Figure15B.BetaPowerBetaOffset,1)), 'LineWidth', 2);
plot(Figure15B.time, squeeze(nanmean(Figure15B.GammaPowerGammaOffset,1)), 'LineWidth', 2);
line([0 0], [-4.5 4.5], 'Color', [0 0 0]); xlim([-1.2 1.2]); ylim([-4.5 4.5]);
xlabel('Time centered on Offset (s)'); ylabel('Power (z-score)'); title('Offset')

subplot(3,5,15);
hold on;
plot(Figure15B.time, squeeze(nanmean(Figure15B.ThetaPowerThetaOnOffset,1)), 'LineWidth', 2);
plot(Figure15B.time, squeeze(nanmean(Figure15B.AlphaPowerAlphaOnOffset,1)), 'LineWidth', 2);
plot(Figure15B.time, squeeze(nanmean(Figure15B.BetaPowerBetaOnOffset,1)), 'LineWidth', 2);
plot(Figure15B.time, squeeze(nanmean(Figure15B.GammaPowerGammaOnOffset,1)), 'LineWidth', 2);
line([0 0], [-4.5 4.5], 'Color', [0 0 0]); xlim([-1.2 1.2]); ylim([-4.5 4.5]);
xlabel('Time centered on On/Offset (s)'); ylabel('Power difference (z-score)'); title('On-Offset')

set(findall(gcf,'-property','FontSize'),'FontSize',18)

legend({'Theta', 'Alpha', 'Beta', 'Gamma'}, 'FontSize', 13, 'location', 'SouthEast'); legend('boxoff');

pn.plotFolder = '/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F_Figures/';
figureName = 'F15B';

saveas(h, [pn.plotFolder, figureName], 'fig');
saveas(h, [pn.plotFolder, figureName], 'epsc');
saveas(h, [pn.plotFolder, figureName], 'png');
