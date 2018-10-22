%% Plot episode characteristics

% pn.data = '/Volumes/EEG/BOSC_Sternberg/B_analyses/B_BOSC/B_eBOSC_180509_full_wl6_noDur/B_180509_MergeRhythmCharacteristics/B_data/A_RhythmCharacteristics/';
% 
% sessions = {'1'; '7'; '8'};
% IDs = {'3130';'3131';'3132';'3134';'3135';'3136';'3138';'3146';'3147';...
%             '3149';'3154';'3156';'3157';'3158';'3159';'3161';'3162';'3233';...
%             '3237';'3239';'3240';'3241';'3242';'3243';'3245';'3248';'3250';...
%             '3251';'3252';'3253';'3255';'3260'};
% 
% for indSession = 1:3
%     for indSubject = 1:numel(IDs)
%         disp(['Session ', num2str(indSession) ' Subject ', num2str(indSubject)])
%         load([pn.data, 'Rhythms_',IDs{indSubject},'_S',sessions{indSession},'_8_15.mat'])
%         RhythmCharacteristicsMerged.CycleDurationMean(indSession,indSubject,:,:)    = RhythmCharacteristics.CycleDurationMean;
%         RhythmCharacteristicsMerged.CycleDurationMax(indSession,indSubject,:,:)     = RhythmCharacteristics.CycleDurationMax;
%         RhythmCharacteristicsMerged.FreqMean(indSession,indSubject,:,:)             = RhythmCharacteristics.FreqMean;
%         RhythmCharacteristicsMerged.FreqMax(indSession,indSubject,:,:)              = RhythmCharacteristics.FreqMax;
%         RhythmCharacteristicsMerged.PowMean(indSession,indSubject,:,:)              = RhythmCharacteristics.PowMean;
%         RhythmCharacteristicsMerged.PowMax(indSession,indSubject,:,:)               = RhythmCharacteristics.PowMax;
%         RhythmCharacteristicsMerged.DetSum(indSession,indSubject,:,:)               = RhythmCharacteristics.DetSum;
%         load([pn.data, 'Bursts_',IDs{indSubject},'_S',sessions{indSession},'_8_15.mat'])
%         BurstCharacteristicsMerged.CycleDurationMean(indSession,indSubject,:,:)    = BurstCharacteristics.CycleDurationMean;
%         BurstCharacteristicsMerged.CycleDurationMax(indSession,indSubject,:,:)     = BurstCharacteristics.CycleDurationMax;
%         BurstCharacteristicsMerged.FreqMean(indSession,indSubject,:,:)             = BurstCharacteristics.FreqMean;
%         BurstCharacteristicsMerged.FreqMax(indSession,indSubject,:,:)              = BurstCharacteristics.FreqMax;
%         BurstCharacteristicsMerged.PowMean(indSession,indSubject,:,:)              = BurstCharacteristics.PowMean;
%         BurstCharacteristicsMerged.PowMax(indSession,indSubject,:,:)               = BurstCharacteristics.PowMax;
%         BurstCharacteristicsMerged.DetSum(indSession,indSubject,:,:)               = BurstCharacteristics.DetSum;
%     end
% end
% 
% %% get timing
% 
% load('/Volumes/EEG/BOSC_Sternberg/B_analyses/B_BOSC/B_eBOSC_180509_full_wl6_noDur/A_180509_full_wl6_noDurThresh/B_data/B_BOSCout/3130_S1_bosc_all.mat')
% 
% Figure14.timevec = bosc_all.time;
% 
% %% sort subjects by rhythmic task SNR
% 
% Figure14.DetSumRhythm = squeeze(nanmedian(nanmedian(RhythmCharacteristicsMerged.DetSum,1),2));
% Figure14.CycleDurationMeanRhythm = squeeze(nanmedian(nanmedian(RhythmCharacteristicsMerged.CycleDurationMean,1),2));
% Figure14.FreqMeanRhythm = squeeze(nanmedian(nanmedian(RhythmCharacteristicsMerged.FreqMean,1),2));
% Figure14.PowMeanRhythm = squeeze(nanmedian(nanmedian(RhythmCharacteristicsMerged.PowMean,1),2));
% 
% Figure14.DetSumBurst = squeeze(nanmedian(nanmedian(BurstCharacteristicsMerged.DetSum,1),2));
% Figure14.CycleDurationMeanBurst = squeeze(nanmedian(nanmedian(BurstCharacteristicsMerged.CycleDurationMean,1),2));
% Figure14.FreqMeanBurst = squeeze(nanmedian(nanmedian(BurstCharacteristicsMerged.FreqMean,1),2));
% Figure14.PowMeanBurst = squeeze(nanmedian(nanmedian(BurstCharacteristicsMerged.PowMean,1),2));
% 
% selectedChans = 44:60; indSession = 3;
% load(['/Volumes/EEG/BOSC_Sternberg/B_analyses/B_BOSC/B_eBOSC_170703_full_wl6/B_extractIndices/B_data/X_8to15_170703_v3.mat']);
% cond = {'L2','L4','L6'}; condFull = {'Load 2','Load 4','Load 6'};
% abn_data = []; eDiffRel = [];
% for c = 1:numel(cond)
%     eDiffRel(:,c) = nanmean(nanmean((X{indSession,1}.([cond{c},'e_amp'])(:,selectedChans,:)-X{indSession,1}.([cond{c},'e_fitBG'])(:,selectedChans,:))./X{indSession,1}.([cond{c},'e_fitBG'])(:,selectedChans,:),3),2)'; 
% end
% [~, SNRsortIDX] = sort(nanmean(eDiffRel,2), 'descend'); % sort by SNR
% 
% Figure14.DetSumRhythm_sorted = squeeze(nanmedian(nanmedian(RhythmCharacteristicsMerged.DetSum(:,SNRsortIDX,44:60,:),1),3));
% Figure14.DetSumBurst_sorted = squeeze(nanmedian(nanmedian(BurstCharacteristicsMerged.DetSum(:,SNRsortIDX,44:60,:),1),3));
% 
% %% save Figure data
% 
% save('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/F14.mat', 'Figure14')

%% load Figure data

addpath('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F4_FigureTools/');

load('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/F14.mat', 'Figure14')

%% Figure: Rhythm vs. burst characteristics

h = figure('units','normalized','position',[.1 .1 .7 .9]); colormap(hot)
subplot(4,2,1);imagesc(Figure14.timevec, [], Figure14.DetSumRhythm, [0 50]); title('Median Detected Episodes');
addTaskTiming(get(gca,'YLim')); xlabel('Time (s)'); ylabel({'Channel';'(anterior-posterior)'}); colorbar; xlim([-2 5])
subplot(4,2,3);imagesc(Figure14.timevec, [], Figure14.CycleDurationMeanRhythm); title('Median Cycle Duration');
addTaskTiming(get(gca,'YLim')); xlabel('Time (s)'); ylabel({'Channel';'(anterior-posterior)'}); colorbar; xlim([-2 5])
subplot(4,2,5);imagesc(Figure14.timevec, [], Figure14.FreqMeanRhythm, [10 12]); title('Median Frequency');
addTaskTiming(get(gca,'YLim')); xlabel('Time (s)'); ylabel({'Channel';'(anterior-posterior)'}); colorbar; xlim([-2 5])
subplot(4,2,7);imagesc(Figure14.timevec, [], Figure14.PowMeanRhythm); title('Median Power');
addTaskTiming(get(gca,'YLim')); xlabel('Time (s)'); ylabel({'Channel';'(anterior-posterior)'}); colorbar; xlim([-2 5])

subplot(4,2,2);imagesc(Figure14.timevec, [], Figure14.DetSumBurst, [0 50]); title('Median Detected Episodes');
addTaskTiming(get(gca,'YLim')); xlabel('Time (s)'); ylabel({'Channel';'(anterior-posterior)'}); colorbar; xlim([-2 5])
subplot(4,2,4);imagesc(Figure14.timevec, [], Figure14.CycleDurationMeanBurst); title('Median Cycle Duration');
addTaskTiming(get(gca,'YLim')); xlabel('Time (s)'); ylabel({'Channel';'(anterior-posterior)'}); colorbar; xlim([-2 5])
subplot(4,2,6);imagesc(Figure14.timevec, [], Figure14.FreqMeanBurst, [10 12]); title('Median Frequency');
addTaskTiming(get(gca,'YLim')); xlabel('Time (s)'); ylabel({'Channel';'(anterior-posterior)'}); colorbar; xlim([-2 5])
subplot(4,2,8);imagesc(Figure14.timevec, [], Figure14.PowMeanBurst); title('Median Power');
addTaskTiming(get(gca,'YLim')); xlabel('Time (s)'); ylabel({'Channel';'(anterior-posterior)'}); colorbar; xlim([-2 5])

set(findall(gcf,'-property','FontSize'),'FontSize',18)
    
pn.plotFolder = ['/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F_Figures/'];
figureName = 'F14A';

saveas(h, [pn.plotFolder, figureName], 'fig');
saveas(h, [pn.plotFolder, figureName], 'epsc');
saveas(h, [pn.plotFolder, figureName], 'png');

%% Figure: plot interindividual differences in cylce duration sustained vs. burst

h = figure('units','normalized','position',[.1 .1 .7 .3]); colormap(hot)
subplot(1,2,1);imagesc(Figure14.timevec, [], Figure14.DetSumRhythm_sorted); title('Rhythmic episodes: Median episode count');
addTaskTiming(get(gca,'YLim')); xlabel('Time (s)'); ylabel('Subjects (sorted by SNR)'); colorbar; xlim([-2 5])

subplot(1,2,2);imagesc(Figure14.timevec, [], Figure14.DetSumBurst_sorted,[0 70]); title('Burst episodes: Median episode count');
addTaskTiming(get(gca,'YLim')); xlabel('Time (s)'); ylabel('Subjects (sorted by SNR)'); colorbar; xlim([-2 5])
set(findall(gcf,'-property','FontSize'),'FontSize',18)

pn.plotFolder = ['/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F_Figures/'];
figureName = 'F14B';

saveas(h, [pn.plotFolder, figureName], 'fig');
saveas(h, [pn.plotFolder, figureName], 'epsc');
saveas(h, [pn.plotFolder, figureName], 'png');
