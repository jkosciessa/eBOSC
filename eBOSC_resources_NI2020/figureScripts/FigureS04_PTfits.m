%% Plot spectra with superimposed power thresholds

% implementation of xaxis labels etc. requires MATLAB R2017+
% original: /Volumes/EEG/BOSC_SternRest/B_scripts_JQK/B_processing_170704/J2_BGoverview_170824.m

% sort subjects by maximum power in eyes closed state, z-score both power
% and thresholds together

% %% some parameters + loading data
% 
% FigureS2.info.channels = 44:60;
% FigureS2.info.freqTask = 2.^[0:.125:6];
% FigureS2.info.freqRest = 2.^[0:.125:6];
% FigureS2.info.sepValue = 7; % seperation between individual traces
% 
% % load background FigureS2.info for resting state
% 
% load(['/Volumes/EEG/BOSC_SternRest/B_analyses/A_eBOSC/A_SternBRest_170703/L_individualThresholds/B_data/Thresholds_170703.mat'], 'Thresholds');
% 
% FigureS2.ThresholdsRest = Thresholds; clear Thresholds;
% 
% % get task spectra
% 
% session = '1';
% pn.bosc = ['/Users/kosciessa/Desktop/mntTardis/Stern_WIP/B_eBOSC_180917_retention/A_eBOSC_180917_retention/B_data/B_BOSCout/eBOSC_S', session, '/'];
% % N = 32 YA
% IDs = {'3130';'3131';'3132';'3134';'3135';'3136';'3138';'3146';'3147';...
%     '3149';'3154';'3156';'3157';'3158';'3159';'3161';'3162';'3233';...
%     '3237';'3239';'3240';'3241';'3242';'3243';'3245';'3248';'3250';...
%     '3251';'3252';'3253';'3255';'3260'};
% fileNames = strcat(IDs, repmat('_bosc.mat',numel(IDs),1));
% IDs = cellfun(@str2num,IDs);
% cond = {'L2','L4','L6'};
% condlab = {'load2', 'load4', 'load6'};
% for indID = 1:length(IDs)
%     % define ID
%     ID = num2str([IDs(indID,1)]); disp(ID);
%     % load BGinfo
%     load([pn.bosc ID, '_bosc.mat'], 'BGinfo');
%     Thresholds.load2_pt(indID,:,:)        = BGinfo.all.load2_pt;
%     Thresholds.load2_pv(indID,:,:)        = BGinfo.all.load2_pv;
%     Thresholds.load2_bg_pow(indID,:,:)    = BGinfo.all.load2_bg_pow;
%     Thresholds.load4_pt(indID,:,:)        = BGinfo.all.load4_pt;
%     Thresholds.load4_pv(indID,:,:)        = BGinfo.all.load4_pv;
%     Thresholds.load4_bg_pow(indID,:,:)    = BGinfo.all.load4_bg_pow;
%     Thresholds.load6_pt(indID,:,:)        = BGinfo.all.load6_pt;
%     Thresholds.load6_pv(indID,:,:)        = BGinfo.all.load6_pv;
%     Thresholds.load6_bg_pow(indID,:,:)    = BGinfo.all.load6_bg_pow;
% end
% 
% FigureS2.ThresholdsTask = Thresholds; clear Thresholds;
% 
% %% save Figure data
% 
% save('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/S2.mat', 'FigureS2')

%% load Figure data

load('/Users/kosciessa/Desktop/eBOSC/figureData/S4.mat', 'FigureS2')

%% plot spectra + superimposed PTs

% sort spectra by descending mean 8-15 Hz power during EC rest
alphaIdx = find(FigureS2.info.freqRest >= 8 & FigureS2.info.freqRest <=15);
IndividualAlphaPowGB_ec = squeeze(mean(mean(FigureS2.ThresholdsRest.ec_bg_pow(:,FigureS2.info.channels,alphaIdx),3),2));
[sortVal, sortInd] = sort(IndividualAlphaPowGB_ec, 'ascend');

h = figure('units','normalized','position',[.1 .1 .7 .7]);

% EYES CLOSED

% z-score both spectra and power thresholds
zData = zscore([squeeze(mean(FigureS2.ThresholdsRest.ec_pt(:,FigureS2.info.channels,:),2)), ...
    squeeze(mean(FigureS2.ThresholdsRest.ec_bg_pow(:,FigureS2.info.channels,:),2))],[],2);

subplot(3,3,[1,4]);
    for count = 1:numel(sortInd)
        indID = sortInd(count);
        hold on; 
        plot(zData(indID,1:size(zData,2)/2)+count*FigureS2.info.sepValue, 'r', 'LineWidth', 2);
        plot(zData(indID,size(zData,2)/2+1:end)+count*FigureS2.info.sepValue, 'k', 'LineWidth', 2);
    end
    xlim([1 numel(FigureS2.info.freqRest)]); ylim([0 count*FigureS2.info.sepValue+FigureS2.info.sepValue]);
    xticks([1:8:numel(FigureS2.info.freqRest)]); xticklabels(round(FigureS2.info.freqRest(xticks),1)); xlabel('Frequency')
    yticks([]); ylabel('Subject')
    title('Eyes Closed');

zData_PT = zData(:,1:49);
zData_BG = zData(:,50:end);
zData_BG_Sub = zData_BG;
zData_BG_Sub(zData_BG>zData_PT) = NaN;
zData_BG_Super = zData_BG;
zData_BG_Super(zData_BG<zData_PT) = NaN;

subplot(3,3,[7]);
hold on;
plot(zData_BG_Sub', 'k', 'LineWidth', 1)
plot(zData_BG_Super', 'r', 'LineWidth', 1)
xlim([1 numel(FigureS2.info.freqRest)]);
xticks([1:8:numel(FigureS2.info.freqRest)]); xticklabels(round(FigureS2.info.freqRest(xticks),1)); xlabel('Frequency')
yticks([]); ylabel('Amplitude (z-score)')
    
% EYES OPEN

% z-score both spectra and power thresholds
zData = zscore([squeeze(mean(FigureS2.ThresholdsRest.eo_pt(:,FigureS2.info.channels,:),2)), ...
    squeeze(mean(FigureS2.ThresholdsRest.eo_bg_pow(:,FigureS2.info.channels,:),2))],[],2);

subplot(3,3,[2,5]);
    for count = 1:numel(sortInd)
        indID = sortInd(count);
        hold on; 
        plot(zData(indID,1:size(zData,2)/2)+count*FigureS2.info.sepValue, 'r', 'LineWidth', 2);
        plot(zData(indID,size(zData,2)/2+1:end)+count*FigureS2.info.sepValue, 'k', 'LineWidth', 2);
    end
    xlim([1 numel(FigureS2.info.freqRest)]); ylim([0 count*FigureS2.info.sepValue+FigureS2.info.sepValue]);
    xticks([1:8:numel(FigureS2.info.freqRest)]); xticklabels(round(FigureS2.info.freqRest(xticks),1)); xlabel('Frequency')
    yticks([]); ylabel('Subject')
    title('Eyes Open');

zData_PT = zData(:,1:49);
zData_BG = zData(:,50:end);
zData_BG_Sub = zData_BG;
zData_BG_Sub(zData_BG>zData_PT) = NaN;
zData_BG_Super = zData_BG;
zData_BG_Super(zData_BG<zData_PT) = NaN;    

subplot(3,3,[8]);
hold on;
plot(zData_BG_Sub', 'k', 'LineWidth', 1)
plot(zData_BG_Super', 'r', 'LineWidth', 1)
xlim([1 numel(FigureS2.info.freqRest)]);
xticks([1:8:numel(FigureS2.info.freqRest)]); xticklabels(round(FigureS2.info.freqRest(xticks),1)); xlabel('Frequency')
yticks([]); ylabel('Amplitude (z-score)')
    
% Task

zData = zscore([squeeze(mean(FigureS2.ThresholdsTask.load6_pt(:,FigureS2.info.channels,:),2)), ...
    squeeze(mean(FigureS2.ThresholdsTask.load6_bg_pow(:,FigureS2.info.channels,:),2))],[],2);

subplot(3,3,[3,6]);
    for count = 1:numel(sortInd)
        indID = sortInd(count);
        hold on; 
        plot(zData(indID,1:size(zData,2)/2)+count*FigureS2.info.sepValue, 'r', 'LineWidth', 2);
        plot(zData(indID,size(zData,2)/2+1:end)+count*FigureS2.info.sepValue, 'k', 'LineWidth', 2);
    end
    xlim([1 numel(FigureS2.info.freqTask)]); ylim([0 count*FigureS2.info.sepValue+FigureS2.info.sepValue]);
    xticks([1:8:numel(FigureS2.info.freqTask)]); xticklabels(round(FigureS2.info.freqTask(xticks),1)); xlabel('Frequency')
    yticks([]); ylabel('Subject')
    set(gca, 'XScale', 'linear')
    title(['Task: Load 6, Session 1']);

zData_PT = zData(:,1:49);
zData_BG = zData(:,50:end);
zData_BG_Sub = zData_BG;
zData_BG_Sub(zData_BG>=zData_PT) = NaN;
zData_BG_Super = zData_BG;
zData_BG_Super(zData_BG<zData_PT) = NaN;    

subplot(3,3,[9]);
hold on;
plot(zData_BG_Sub', 'k', 'LineWidth', 1)
plot(zData_BG_Super', 'r', 'LineWidth', 1)
xlim([1 numel(FigureS2.info.freqRest)]);
xticks([1:8:numel(FigureS2.info.freqRest)]); xticklabels(round(FigureS2.info.freqRest(xticks),1)); xlabel('Frequency')
yticks([]); ylabel('Amplitude (z-score)')
    
suptitle('Spectra and power thresholds');
set(findall(gcf,'-property','FontSize'),'FontSize',18)

pn.plotFolder = '/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F_Figures/';
figureName = 'S2';

saveas(h, [pn.plotFolder, figureName], 'fig');
saveas(h, [pn.plotFolder, figureName], 'epsc');
saveas(h, [pn.plotFolder, figureName], 'png');
