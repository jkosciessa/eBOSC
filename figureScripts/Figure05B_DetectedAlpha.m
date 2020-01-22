% IMPORTANT: to plot trial-wise power within and without Figure6B.detected segments,
% run B_BOSC/J_temporalOscillationDistribution_400ms_170518 first!

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
% %% load trial-wise indication matrix of when rhythms were Figure6B.detected
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
% Figure6B.timeVec = timeBOSC(idx.BOSCbegin:idx.BOSCend);
% 
% %% get sorting of maxima
% 
% alphaChannel = 60; % plot from 02
% 
% % % get maximum for each subject
% % for indSub = 1:32
% %     % concatenate data from the three conditions
% %     tmp_data = cat(1,trialAlpha.load2{indSub, alphaChannel}, trialAlpha.load4{indSub, alphaChannel}, trialAlpha.load6{indSub, alphaChannel});
% %     tmp_data = squeeze(tmp_data(:,idx.begin+1:idx.end));
% %     mxData = max(max(tmp_data));
% %     colorMax(1,indSub) = mxData;
% %     clear tmp_data;
% % end
% % [~, sortInd] = sort(colorMax, 'descend');
% % 
% % sortSelect = sortInd([1,5,9,23,27,30,32]);
% 
% %% prepare data for representative subjects
% 
% selectedSubjects = [4,5,32,24,9,22];
% 
% Figure6B.detected = []; Figure6B.undetected = []; Figure6B.mxData = [];
% for indSubject = 1:numel(selectedSubjects)
%     indSub = selectedSubjects(indSubject);
%     % concatenate data from the three conditions
%     TrialData = cat(1,trialAlpha.load2{indSub, alphaChannel}, trialAlpha.load4{indSub, alphaChannel}, trialAlpha.load6{indSub, alphaChannel});
%     %TrialData = cat(1,trialAlpha.load2{indSub, alphaChannel});
%     TrialData = squeeze(TrialData(:,idx.begin:idx.end));
%     % z-standardize each trial
%     for indTrial = 1:size(TrialData,1)
%         TrialData(indTrial,:) = zscore(TrialData(indTrial,:));
%     end
%     tmp_detected = squeeze(AlphaTrials(indSub,1:size(TrialData,1),idx.BOSCbegin:idx.BOSCend));
%     Figure6B.mxData{indSubject} = max(max(TrialData));
%     % Figure6B.detected segments on the left; Figure6B.undetected time points on the right
%     Figure6B.detected{indSubject} = tmp_detected(:,1:size(TrialData,2)).*TrialData; 
%     Figure6B.detected{indSubject}(Figure6B.detected{indSubject}==0) = NaN;
%     Figure6B.undetected{indSubject} = (1-tmp_detected(:,1:size(TrialData,2))).*TrialData; 
%     Figure6B.undetected{indSubject}(Figure6B.undetected{indSubject}==0) = NaN;
% end
% 
% save('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/F6B.mat', 'Figure6B')

%% plot detected vs. non-detected power by trial

load('/Users/kosciessa/Desktop/eBOSC/figureData/F5B.mat', 'Figure6B')

hfig = figure('units','normalized','position',[.1 .1 .7 .7]);
[ha, pos] = tight_subplot(7,1,[.01 .01],[.005 .005],[.005 .005]);
for indPlot = 1:6
    axes(ha(indPlot)); 
    [hax1, cax] = imagescwithnan([Figure6B.detected{indPlot}, Figure6B.undetected{indPlot}], hot, [0 0 0]);
    line([numel(Figure6B.timeVec) numel(Figure6B.timeVec)],get(gca,'XLim'), ...
        'LineStyle', '-', 'Color', 'w', 'LineWidth', 3)
    colorbar('off')
    ax{indPlot} = gca;
    ax{indPlot}.LineWidth = 10;
    colorMax(1,indPlot) = Figure6B.mxData{indPlot};
    ax{indPlot}.XTick = [];
    ax{indPlot}.YTick = [];
    clear tmp_data;
end
% add schematic of time periods
ha(7).Position = [ha(7).Position(1),ha(7).Position(2)+.07,ha(7).Position(3),ha(7).Position(4)-.07];
axes(ha(7));
timeComplete = [Figure6B.timeVec,5+2.3960+Figure6B.timeVec];
plot(timeComplete, NaN(1,numel(timeComplete)));
patches.timeVec = [-2.2 -1.4 -1.2 -.2 0 3 3.2, 5+2.3960+[-1.4 -1.2 -.2 0 3 3.2]];
patches.timeVecActual = [-2.2 -1.4 -1.2 -.2 0 3 3.2, [-1.4 -1.2 -.2 0 3 3.2]];
patches.colorVec = [1 1 1;0 0 0; 1 1 1; 0 0 0; 1 0 0; 0 0 0; [1 1 1;0 0 0; 1 1 1;0 0 0; 1 0 0;0 0 0]];
for indP = 1:numel(patches.timeVec)-1
    p = patch([patches.timeVec(indP) patches.timeVec(indP+1) patches.timeVec(indP+1) patches.timeVec(indP)], ...
                [.08 .08 [.90 .90]], patches.colorVec(indP,:));
    p.EdgeColor = 'none';
end; hold on;
% add individual yellow borders
xlim([min(Figure6B.timeVec) max(timeComplete)]);
colorbar('off')
ax{7} = gca;
ax{7}.LineWidth = 10;
ax{7}.XTick = [-1.4 -1.2 -.2 0 3 3.2];
ax{7}.YTick = [];
line([5,5],[0 1], 'LineStyle', '-', 'Color', 'k', 'LineWidth', 3)
for indSub = 1:6
    tmp_color = sort(colorMax, 'descend');
    tmp_color = tmp_color(end-5:end);
    tmp_color = tmp_color./max(tmp_color);
    curCol = tmp_color(indSub);
    ax{indSub}.XColor = [curCol curCol 0];
    ax{indSub}.YColor = [curCol curCol 0];
    ax{indSub}.ZColor = [curCol curCol 0];
end
ax{7}.XColor = [0 0 0];
ax{7}.YColor = [0 0 0];
ax{7}.ZColor = [0 0 0];
set(gca,'TickLength',[0 0])
set(gca, 'XTick', patches.timeVec([3,5,6,9,11,12]))
set(gca, 'XTickLabel', patches.timeVecActual([3,5,6,9,11,12]))
set(gca, 'XAxisLocation', 'bottom')
set(gca, 'TickDir', 'in')
xlabel('Time relative to Retention onset (seconds)')
set(findall(gcf,'-property','FontSize'),'FontSize',20)

pn.plotFolder = '/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F_Figures/';

saveas(hfig, [pn.plotFolder, 'F2B_v2'], 'epsc');
saveas(hfig, [pn.plotFolder, 'F2B_v2'], 'fig');
saveas(hfig, [pn.plotFolder, 'F2B_v2'], 'png');
%close(hfig);

%% plot lower part of Figure 6

load('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/F6A.mat', 'Figure6A')
addpath('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F4_FigureTools')

restoredefaultpath
pn.shadedError = ['/Volumes/EEG/BOSC_Sternberg/T_tools/shadedErrorBar-7002ebc']; addpath(pn.shadedError);

hfig = figure('units','normalized','position',[.1 .1 .7 .2]);
subplot(1,2,1); 
    patches.timeVec = [-1.4, -1.2; -.2 0; 3, 3.2];
    patches.colorVec = [.8 .8 .8; .8 .8 .8; .6 .6 .6];
    for indP = 1:size(patches.timeVec,1)
        YLim = [-.25 .5];
        p = patch([patches.timeVec(indP,1) patches.timeVec(indP,2) patches.timeVec(indP,2) patches.timeVec(indP,1)], ...
                    [YLim(1) YLim(1)  YLim(2), YLim(2)], patches.colorVec(indP,:));
        p.EdgeColor = 'none';
    end; hold on;
    yyaxis left; set(gca, 'ylim', [-.25 .5]); %plot(Figure6A.timeVec,nanmean(Figure6A.SubAvg.detTF_z,1), 'k', 'LineWidth', 2); %hold on; plot(Figure6A.timeVec,nanmean(SubAvg.z_Power,1), '--k'); 
    standError = nanstd(Figure6A.SubAvg.detTF_z,1)./sqrt(size(Figure6A.SubAvg.detTF_z,1));
    l1 = shadedErrorBar(Figure6A.timeVec,nanmean(Figure6A.SubAvg.detTF_z,1),standError, 'lineprops', {'color', [0 0 0],'linewidth', 3}, 'patchSaturation', .05);
    ylabel({'rhythmic power';'(z-score)'},'FontSize',15);
    yyaxis right; set(gca, 'ylim', [0 40]); %plot(Figure6A.timeVec,nanmean(Figure6A.SubAvg.det,1)*100,'r', 'LineWidth', 2); 
    ylabel('% trials rhythmic','FontSize',15, 'Color', 'r'); set(gca, 'YTickLabel', [0 10 20 30 40], 'YColor','r'); xlabel('Time (s)');xlim([Figure6A.timeVec(1),Figure6A.timeVec(end)]);
    standError = nanstd(Figure6A.SubAvg.det.*100,1)./sqrt(size(Figure6A.SubAvg.det,1));
    l2 = shadedErrorBar(Figure6A.timeVec,nanmean(Figure6A.SubAvg.det.*100,1),standError, 'lineprops', {'color', [1 0 0],'linewidth', 3}, 'patchSaturation', .05);
    xlim([-2.4 5])
    text(.8,10, 'Retention')
subplot(1,2,2); cla;
    for indP = 1:size(patches.timeVec,1)
        YLim = [-.8 .35];
        p = patch([patches.timeVec(indP,1) patches.timeVec(indP,2) patches.timeVec(indP,2) patches.timeVec(indP,1)], ...
                    [YLim(1) YLim(1) YLim(2), YLim(2)], patches.colorVec(indP,:));
        p.EdgeColor = 'none';
    end; hold on;
    % add error bars
    yyaxis left; set(gca, 'ylim', [-.8 .35]);
    standError = nanstd(Figure6A.SubAvg.undetTF_z,1)./sqrt(size(Figure6A.SubAvg.undetTF_z,1));
    l3 = shadedErrorBar(Figure6A.timeVec,nanmean(Figure6A.SubAvg.undetTF_z,1),standError, 'lineprops', {'color', [0 0 0],'linewidth', 3}, 'patchSaturation', .05);
    ylabel({'arrhythmic power';'(z-score)'},'FontSize',15);
    yyaxis right; set(gca, 'ylim', [0 40]); %plot(Figure6A.timeVec,nanmean(Figure6A.SubAvg.det,1)*100,'r', 'LineWidth', 2); 
    ylabel('% trials rhythmic','FontSize',15, 'Color', 'r'); set(gca, 'YTickLabel', [0 10 20 30 40], 'YColor','r'); xlabel('Time (s)'); xlim([Figure6A.timeVec(1),Figure6A.timeVec(end)]);
    standError = nanstd(Figure6A.SubAvg.det.*100,1)./sqrt(size(Figure6A.SubAvg.det,1));
    l4 = shadedErrorBar(Figure6A.timeVec,nanmean(Figure6A.SubAvg.det.*100,1),standError, 'lineprops', {'color', [1 0 0],'linewidth', 3}, 'patchSaturation', .05);
    xlim([-2.4 5])
    text(.8,10, 'Retention')
set(findall(gcf,'-property','FontSize'),'FontSize',22)

pn.plotFolder = '/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F_Figures/';

saveas(hfig, [pn.plotFolder, 'F6B_lower'], 'epsc');
saveas(hfig, [pn.plotFolder, 'F6B_lower'], 'fig');
saveas(hfig, [pn.plotFolder, 'F6B_lower'], 'png');