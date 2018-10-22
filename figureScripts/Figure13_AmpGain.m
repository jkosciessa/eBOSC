% - used to be /Volumes/EEG/BOSC_SternRest/B_scripts_JQK/B_processing_170704/E_JQKplotAmpAbn_170816B.mlx

% pn.root = '/Volumes/EEG/BOSC_Sternberg/B_analyses/B_BOSC/B_eBOSC_170703_full_wl6/';
% pn.tools = [pn.root, 'C_interrelations/T_tools/']; addpath(genpath(pn.tools));
% pn.XTask = [pn.root, 'B_extractIndices/B_data/'];
% pn.XRest = ['/Volumes/EEG/BOSC_SternRest/B_analyses/A_eBOSC/A_SternBRest_170703/B_extractIndices/B_data/'];
% pn.plotFolder = '/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F_Figures/';
% 
% Figure13.colorm = brewermap(4,'RdYlBu'); % use external function to create colormap
% 
% info.channels = 44:60;
% 
% load([pn.XRest, 'X15B_8to15_170703.mat']);
% 
% cond = {'EC1','EC2','EO1','EO2'};
% 
% for c = 1:numel(cond)
%     Figure13.Rest.abn_data(:,c) = nanmean(nanmean(X{1,1}.([cond{c},'e_abn'])(:,info.channels,:),3),2)';
%     % xAmp-BG
%     Figure13.Rest.oDiff(:,c) = nanmean(nanmean(X{1,1}.([cond{c},'o_amp'])(:,info.channels,:) - X{1,1}.([cond{c},'o_fitBG'])(:,info.channels,:),3),2)';
%     Figure13.Rest.eDiff(:,c) = nanmean(nanmean(X{1,1}.([cond{c},'e_amp'])(:,info.channels,:) - X{1,1}.([cond{c},'e_fitBG'])(:,info.channels,:),3),2)';
% end
% 
% %% get task data
% 
% pn.root = '/Volumes/EEG/BOSC_Sternberg/B_analyses/B_BOSC/B_eBOSC_170703_full_wl6/';
% pn.tools = [pn.root, 'C_interrelations/T_tools/']; addpath(genpath(pn.tools));
% pn.Xdata = [pn.root, 'B_extractIndices/B_data/'];
% 
% sessions = {'1'; '7'; '8'};
% indSession = 1;
% info.channels = 44:60;
% cond = {'L2','L4','L6'};
% 
% load([pn.XTask, 'X_8to15_170703_v3.mat']);
% 
% for c = 1:numel(cond)
%     Figure13.Task.abn_data(:,c) = nanmean(nanmean(X{1,1}.([cond{c},'e_abn'])(:,info.channels,:),3),2)';
%     % xAmp-BG
%     Figure13.Task.oDiff(:,c) = nanmean(nanmean(X{1,1}.([cond{c},'o_amp'])(:,info.channels,:) - X{1,1}.([cond{c},'o_fitBG'])(:,info.channels,:),3),2)';
%     Figure13.Task.aDiff(:,c) = nanmean(nanmean(X{1,1}.([cond{c},'a_amp'])(:,info.channels,:) - X{1,1}.([cond{c},'a_fitBG'])(:,info.channels,:),3),2)';
%     Figure13.Task.eDiff(:,c) = nanmean(nanmean(X{1,1}.([cond{c},'e_amp'])(:,info.channels,:) - X{1,1}.([cond{c},'e_fitBG'])(:,info.channels,:),3),2)';
% end
% 
% %% save Figure data
% 
% save('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/F13.mat', 'Figure13')

%% load Figure data

load('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/F13.mat', 'Figure13')

%% plot arrhythmic bias in both rest and task
% ((eAmp-eBG)-(oAmp-oBG))./(eAmp-eBG)

h = figure('units','normalized','position',[.1 .1 .7 .4]);
rho = []; p = [];
[rho p] = corr(1-Figure13.Rest.abn_data(:,:),(Figure13.Rest.eDiff(:,:)-...
    Figure13.Rest.oDiff(:,:))./Figure13.Rest.eDiff(:,:), 'type', 'Pearson', 'rows', 'complete')
scientificVals = cell(size(p));
for indP = 1:numel(p)
    n = p(indP);
    n_exp=floor(log10(n));
    n_coef=n/(abs(1*10^n_exp));
    coef = sprintf('%.0f', round(n_coef,2));
    exp = sprintf('%.0f', round(n_exp,0));
    if abs(n_exp)<= 3
        scientificVals{indP} = sprintf('%.3f', round(p(indP),3));
    else
        scientificVals{indP} = [coef, 'e^', exp];
    end
end

subplot(1,2,1);
for indCond = 1:4
    scatter(1-Figure13.Rest.abn_data(:,indCond), (Figure13.Rest.eDiff(:,indCond)-...
        Figure13.Rest.oDiff(:,indCond))./Figure13.Rest.eDiff(:,indCond),[],Figure13.colorm(indCond,:),'filled');
    hold on;
end
hold on; plot([1 0],[1 0],'k')
axis([0 1 0 1])
xlabel('Arrhythmic duration (1-abundance)'); ylabel('Detection-induced rhythmic amplitude gain')
title('Resting data')
leg{1} = legend(['EC: r = ', num2str(rho(1,1),'%.2f'), '; p = ', scientificVals{1,1}],...
    ['ECEO - EC: r = ', num2str(rho(2,2),'%.2f'), '; p = ', scientificVals{2,2}]',...
    ['EO: r = ', num2str(rho(3,3),'%.2f'), '; p = ', scientificVals{3,3}]',...
    ['ECEO - EO: r = ', num2str(rho(4,4),'%.2f'), '; p = ', scientificVals{4,4}]',...
    'location','northwest', 'Orientation', 'vertical');
legend('boxoff');

%% plot bias in task

% ((eAmp-eBG)-(oAmp-oBG))./(eAmp-eBG)

[rho p] = corr(1-Figure13.Task.abn_data(:,:),(Figure13.Task.eDiff(:,:)-...
    Figure13.Task.oDiff(:,:))./Figure13.Task.eDiff(:,:), 'type', 'Pearson', 'rows', 'complete')
scientificVals = cell(size(p));
for indP = 1:numel(p)
    n = p(indP);
    n_exp=floor(log10(n));
    n_coef=n/(abs(1*10^n_exp));
    coef = sprintf('%.0f', round(n_coef,2));
    exp = sprintf('%.0f', round(n_exp,0));
    if abs(n_exp)<= 3
        scientificVals{indP} = sprintf('%.3f', round(p(indP),3));
    else
        scientificVals{indP} = [coef, 'e^', exp];
    end
end

cond = {'L2','L4','L6'};

subplot(1,2,2);
for indCond = 1:numel(cond)
    scatter(1-Figure13.Task.abn_data(:,indCond), (Figure13.Task.eDiff(:,indCond)-...
        Figure13.Task.oDiff(:,indCond))./Figure13.Task.eDiff(:,indCond),[],Figure13.colorm(indCond,:),'filled');
    hold on;
end
hold on; plot([1 0],[1 0],'k')
axis([0 1 0 1])
xlabel('Arrhythmic duration (1-abundance)'); ylabel('Detection-induced rhythmic amplitude gain')
title('Task data: Retention period')
leg{2} = legend(['Load 2: r = ', num2str(rho(1,1),'%.2f'), '; p = ', scientificVals{1,1}],...
    ['Load 4: r = ', num2str(rho(2,2),'%.2f'), '; p = ', scientificVals{2,2}]',...
    ['Load 6: r = ', num2str(rho(3,3),'%.2f'), '; p = ', scientificVals{3,3}]',...
    'location','northwest', 'Orientation', 'vertical');
legend('boxoff');

set(findall(gcf,'-property','FontSize'),'FontSize',16)

h_sup = suptitle('Arrhythmic content linearly decreases rhythmic amplitude estimates');
set(h_sup, 'FontWeight','bold', 'FontSize', 18);

set(leg{1}, 'FontSize', 12);
set(leg{2}, 'FontSize', 12);

%% save plot

saveas(h, [pn.plotFolder, 'F13'], 'fig');
saveas(h, [pn.plotFolder, 'F13'], 'epsc');
saveas(h, [pn.plotFolder, 'F13'], 'png');
