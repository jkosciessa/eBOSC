% %% Plot time-wise detected vs. undetected power
% 
% clear all; clc;
% 
% sessions = {'1', '7', '8'};
% indSession = 1;
% 
% pn.root = '/Volumes/EEG/BOSC_Sternberg/B_analyses/B_BOSC/B_eBOSC_170703_full_wl6/J_AbnPowDynamics/A_AlphaTheta/';
% 
% addpath(genpath('/Volumes/EEG/BOSC_Sternberg/T_tools/tight_subplot'))
% addpath('/Volumes/EEG/BOSC_Sternberg/T_tools/fnct_JQK/');
% 
% %% load time-frequency data
% 
% load([pn.root, 'B_data/J3_trialAlphaTheta_S',sessions{indSession},'_170703.mat'], 'trialAlpha', 'time');
% 
% %% load trial-wise indication matrix of when rhythms were Figure6A.detected
% 
% load([pn.root, 'B_data/J2_AlphaThetaDetected_S',sessions{indSession},'_170703.mat'], 'AlphaTrials')
% 
% %% align BOSC timing with TFR timing
% 
% idx.begin = find(time > -2.4, 1, 'first');
% idx.end = find(time <= 5, 1, 'last'); % note: only goes until 6 seconds
% 
% numSamples = idx.end-idx.begin;
% 
% timeBOSC = -4:0.004:7; timeBOSC(1) = [];
% 
% idx.BOSCbegin = find(timeBOSC > -2.4, 1, 'first');
% idx.BOSCend = find(timeBOSC <= 5, 1, 'last');
% 
% Figure6A.timeVec = timeBOSC(idx.BOSCbegin:idx.BOSCend);
% 
% %% alpha load 6 version with scaling
% 
% alphaChannel = 44:60;
% 
% % get maximum for each subject
% for indSub = 1:32
%     DataChanAvg = cat(3, trialAlpha.load6{indSub, 44:60});
%     tmp_data = nanmean(DataChanAvg,3); clear DataChanAvg;
%     tmp_data = squeeze(tmp_data(:,idx.begin+1:idx.end));
%     mxData = max(max(tmp_data));
%     colorMax(1,indSub) = mxData;
%     clear tmp_data;
% end
% [~, sortInd] = sort(colorMax, 'descend');
% 
% replaceWith = 0;
% 
% indRow = 1; SubAvg = [];
% for indSub = sortInd
%     TrialData_z = []; TrialData = [];
%     % average across posterior channels
%     DataChanAvg = cat(3, trialAlpha.load6{indSub, 44:60});
%     TrialData = nanmean(DataChanAvg,3); clear DataChanAvg;
%     load6begin = size(trialAlpha.load2{indSub, 1},1)+size(trialAlpha.load4{indSub, 1},1)+1;
%     TrialData = squeeze(TrialData(:,idx.begin:idx.end));
%     % z-standardize each trial
%     for indTrial = 1:size(TrialData,1)
%         TrialData_z(indTrial,:) = zscore(TrialData(indTrial,:));
%     end
%     tmp_detected = squeeze(AlphaTrials(indSub,load6begin:load6begin+size(TrialData,1)-1,idx.BOSCbegin:idx.BOSCend));
%     % detected segments on the left; undetected time points on the right
%     RawPower = TrialData;
%     z_Power = TrialData_z;
%     detected_z = tmp_detected(:,1:size(TrialData_z,2)).*TrialData_z; detected_z(detected_z==0) = replaceWith;
%     undetected_z = (1-tmp_detected(:,1:size(TrialData_z,2))).*TrialData_z; undetected_z(undetected_z==0) = replaceWith;
%     Figure6A.SubAvg.RawPower(indRow,:) = nanmean(RawPower,1);
%     Figure6A.SubAvg.z_Power(indRow,:) = nanmean(z_Power,1);
%     Figure6A.SubAvg.det(indRow,:) = nanmean(tmp_detected,1);
%     Figure6A.SubAvg.undet(indRow,:) = nanmean(1-tmp_detected,1);
%     Figure6A.SubAvg.detTF_z(indRow,:) = nanmean(detected_z,1);
%     Figure6A.SubAvg.undetTF_z(indRow,:) = nanmean(undetected_z,1);
%     indRow = indRow + 1;
% end
% 
% %% save Figure data
% 
% save('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/F6A.mat', 'Figure6A')

%% load Figure data

load('/Users/kosciessa/Desktop/eBOSC/figureData/F5A.mat', 'Figure6A')
addpath('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F4_FigureTools')

%% plot Figure

h = figure('units','normalized','position',[.1 .1 .7 .6]);
subplot(3,2,1); imagesc(Figure6A.timeVec,1:32,zscore(Figure6A.SubAvg.RawPower,[],2)); title('Z-Stnd. TF Power','FontSize',15); ylabel('Subject'); xlabel('Time');xlim([Figure6A.timeVec(1),Figure6A.timeVec(end)]); c = colorbar; ylabel(c, 'Power [z-score]');
addTaskTiming(get(gca,'YLim'))
subplot(3,2,2); imagesc(Figure6A.timeVec,1:32,Figure6A.SubAvg.det*100, [0 100]); title('% rhythmic trials','FontSize',15); ylabel('Subject'); xlabel('Time');xlim([Figure6A.timeVec(1),Figure6A.timeVec(end)]); c = colorbar; ylabel(c, 'Rhythmic trials [%]');
addTaskTiming(get(gca,'YLim'))
subplot(3,2,3); imagesc(Figure6A.timeVec,1:32,Figure6A.SubAvg.detTF_z, [-.5 1]); title('z-Stnd. TF Power: Detected','FontSize',15); ylabel('Subject'); xlabel('Time');xlim([Figure6A.timeVec(1),Figure6A.timeVec(end)]); c = colorbar; ylabel(c, 'Power [z-score]');
addTaskTiming(get(gca,'YLim'))
subplot(3,2,4); imagesc(Figure6A.timeVec,1:32,Figure6A.SubAvg.undetTF_z, [-.5 1]); title('z-Stnd. TF Power: Undetected','FontSize',15); ylabel('Subject'); xlabel('Time');xlim([Figure6A.timeVec(1),Figure6A.timeVec(end)]); c = colorbar; ylabel(c, 'Power [z-score]');
addTaskTiming(get(gca,'YLim'))
subplot(3,2,5); 
patches.timeVec = [-1.4, -1.2; -.2 0; 3, 3.2];
patches.colorVec = [.8 .8 .8; .8 .8 .8; .6 .6 .6];
for indP = 1:size(patches.timeVec,1)
    YLim = [-.25 .5];
    p = patch([patches.timeVec(indP,1) patches.timeVec(indP,2) patches.timeVec(indP,2) patches.timeVec(indP,1)], ...
                [YLim(1) YLim(1)  YLim(2), YLim(2)], patches.colorVec(indP,:));
    p.EdgeColor = 'none';
end; hold on;
yyaxis left; set(gca, 'ylim', [-.25 .5]); plot(Figure6A.timeVec,nanmean(Figure6A.SubAvg.detTF_z,1), 'k', 'LineWidth', 2); %hold on; plot(Figure6A.timeVec,nanmean(SubAvg.z_Power,1), '--k'); 
ylabel('z-Stnd. TF Power rhythmic','FontSize',15);
yyaxis right; set(gca, 'ylim', [0 40]); plot(Figure6A.timeVec,nanmean(Figure6A.SubAvg.det,1)*100,'r', 'LineWidth', 2); ylabel('% trials rhythmic','FontSize',15, 'Color', 'r'); set(gca, 'YTickLabel', [0 10 20 30 40], 'YColor','r'); xlabel('Time');xlim([Figure6A.timeVec(1),Figure6A.timeVec(end)]);
subplot(3,2,6); 
for indP = 1:size(patches.timeVec,1)
    YLim = [-.8 .35];
    p = patch([patches.timeVec(indP,1) patches.timeVec(indP,2) patches.timeVec(indP,2) patches.timeVec(indP,1)], ...
                [YLim(1) YLim(1) YLim(2), YLim(2)], patches.colorVec(indP,:));
    p.EdgeColor = 'none';
end; hold on;
yyaxis left; set(gca, 'ylim', [-.8 .35]);plot(Figure6A.timeVec,nanmean(Figure6A.SubAvg.undetTF_z,1), 'LineWidth', 2);%hold on; plot(Figure6A.timeVec,nanmean(SubAvg.z_Power,1), '--k'); 
ylabel('z-Stnd. TF Power arrhythmic','FontSize',15);
yyaxis right; set(gca, 'ylim', [0 40]); plot(Figure6A.timeVec,nanmean(Figure6A.SubAvg.det,1)*100,'r', 'LineWidth', 2);  ylabel('% trials rhythmic','FontSize',15, 'Color', 'r'); set(gca, 'YTickLabel', [0 10 20 30 40], 'YColor','r'); xlabel('Time'); xlim([Figure6A.timeVec(1),Figure6A.timeVec(end)]);

set(findall(gcf,'-property','FontSize'),'FontSize',15)

pn.plotFolder = '/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F_Figures/';

colormap(hot);

saveas(h, [pn.plotFolder, 'F6A'], 'epsc');
saveas(h, [pn.plotFolder, 'F6A'], 'fig');
saveas(h, [pn.plotFolder, 'F6A'], 'png');

%% task timing function

function addTaskTiming(limits)
    hold on;
    line([3 3],limits,'LineWidth', 1, 'Color',[1 1 1], 'LineStyle', '-'); 
    line([3.2 3.2],limits,'LineWidth', 1, 'Color',[1 1 1], 'LineStyle', '-'); 
    line([0 0],limits,'LineWidth', 1, 'Color',[1 1 1], 'LineStyle', '-'); 
    line([-.2 -.2],limits,'LineWidth', 1, 'Color',[1 1 1], 'LineStyle', '-')
    line([-1.2 -1.2],limits,'LineWidth', 1, 'Color',[1 1 1], 'LineStyle', '-'); 
    line([-1.4 -1.4],limits,'LineWidth', 1, 'Color',[1 1 1], 'LineStyle', '-')
end