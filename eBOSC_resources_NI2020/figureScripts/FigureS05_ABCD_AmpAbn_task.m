%% Plot amplitude-abundance relationships dueing task

% - used to be /Volumes/FB-LIP/ConMem/EEG/BOSC_SternRest/B_scripts_JQK/B_processing_170704/E_JQKplotAmpAbn_170816B_Fig8.m

% %% get data
% 
% pn.root = '/Volumes/EEG/BOSC_Sternberg/B_analyses/B_BOSC/B_eBOSC_170703_full_wl6/';
% pn.tools = [pn.root, 'C_interrelations/T_tools/']; addpath(genpath(pn.tools));
% pn.Xdata = [pn.root, 'B_extractIndices/B_data/'];
% 
% %% setup
% 
% sessions = {'1'; '7'; '8'};
% indSession = 1;
% info.channels = 44:60;
% 
% addpath('/Volumes/EEG/BOSC_SternRest/B_analyses/A_eBOSC/A_SternBRest_170703/C_interrelations/T_tools/brewermap/')
% FigureS3.colorm = brewermap(4,'RdYlBu'); % use external function to create colormap
% 
% %% load simulated data
% 
% FigureS3.overallSim = load(['/Volumes/EEG/BOSC_Simulation/17_JQK_RED/D_data/D_simulation/_SimDataforRestingStateComparison_overall_v4.mat'], 'simData');
% 
% %% load task data
% 
% load([pn.Xdata, 'X_8to15_170703_v3.mat']);
% 
% FigureS3.cond = {'L2','L4','L6'};
% 
% for c = 1:numel(FigureS3.cond)
%     FigureS3.e_BGfit(:,c) = nanmean(nanmean(X{indSession,1}.([FigureS3.cond{c},'e_fitBG'])(:,info.channels,:),3),2)';
%     FigureS3.abn_data(:,c) = nanmean(nanmean(X{indSession,1}.([FigureS3.cond{c},'e_abn'])(:,info.channels,:),3),2)';
%     % xAmp./BG
%     FigureS3.oRel(:,c) = nanmean(nanmean(X{indSession,1}.([FigureS3.cond{c},'o_amp'])(:,info.channels,:)./X{indSession,1}.([FigureS3.cond{c},'o_fitBG'])(:,info.channels,:),3),2)';
%     FigureS3.eDiffRel(:,c) = nanmean(nanmean((X{indSession,1}.([FigureS3.cond{c},'e_amp'])(:,info.channels,:)-X{indSession,1}.([FigureS3.cond{c},'e_fitBG'])(:,info.channels,:))./X{indSession,1}.([FigureS3.cond{c},'e_fitBG'])(:,info.channels,:),3),2)';
% end
% 
% %% save Figure data
% 
% save('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/S3.mat', 'FigureS3')

%% load Figure data

load('/Users/kosciessa/Desktop/eBOSC/figureData/S5.mat', 'FigureS3')

%% Plot

h = figure('units','normalized','position',[.1 .1 .8 .3], 'DefaultAxesFontSize',15);

%% SNR (overall, NOT rhythm-exclusive)

% Calculate separate correlations fits below an empirical SNR of 6 and above (for all data here).

% below 6

tmp_abn = reshape(FigureS3.abn_data, [],1);
tmp_snr = reshape(FigureS3.oRel, [],1);

[idx, ~] = find(tmp_snr < 5);

tmp_abn_bel = tmp_abn(idx);
tmp_snr_bel = tmp_snr(idx);

[rho_bel p_bel] = corr(tmp_abn_bel, tmp_snr_bel, 'type', 'Pearson', 'rows', 'complete')
scientificVals_bel = cell(size(p_bel));
for indP = 1:numel(p_bel)
    n = p_bel(indP);
    n_exp=floor(log10(n));
    n_coef=n/(abs(1*10^n_exp));
    coef = sprintf('%.0f', round(n_coef,2));
    exp = sprintf('%.0f', round(n_exp,0));
    if abs(n_exp)<= 3
        scientificVals_bel{indP} = sprintf('%.3f', round(p_bel(indP),3));
    else
        scientificVals_bel{indP} = [coef, 'e^', exp];
    end
end

[idx, ~] = find(tmp_snr >= 5);

tmp_abn_abv = tmp_abn(idx);
tmp_snr_abv = tmp_snr(idx);

[rho_abv p_abv] = corr(tmp_abn_abv, tmp_snr_abv, 'type', 'Pearson', 'rows', 'complete')
scientificVals_abv = cell(size(p_abv));
for indP = 1:numel(p_abv)
    n = p_abv(indP);
    n_exp=floor(log10(n));
    n_coef=n/(abs(1*10^n_exp));
    coef = sprintf('%.0f', round(n_coef,2));
    exp = sprintf('%.0f', round(n_exp,0));
    if abs(n_exp)<= 3
        scientificVals_abv{indP} = sprintf('%.3f', round(p_abv(indP),3));
    else
        scientificVals_abv{indP} = [coef, 'e^', exp];
    end
end

subplot(1,3,1)
for indCond = 1:numel(FigureS3.cond)
    scatter(FigureS3.abn_data(:,indCond),(FigureS3.oRel(:,indCond)),[],FigureS3.colorm(indCond,:),'filled');
    hold on;
end
hold on; plot(FigureS3.overallSim.simData.Abn, FigureS3.overallSim.simData.SNR, 'k', 'lineWidth', 1.5);
Fit = polyfit(tmp_abn_bel,tmp_snr_bel,1);
h1 = plot(tmp_abn_bel, polyval(Fit,tmp_abn_bel));
Fit = polyfit(tmp_abn_abv,tmp_snr_abv,1);
h2 = plot(tmp_abn_abv, polyval(Fit,tmp_abn_abv));
xlabel('Abundance'); ylabel('SNR [Amplitude/background]')
h_legend1 = legend([h2, h1], ['SNR above 6: r = ', num2str(rho_abv,'%.2f'), '; p = ', scientificVals_abv{1,1}], ...
    ['SNR below 6: r = ', num2str(rho_bel,'%.2f'), '; p = ', scientificVals_bel{1,1}]', ...
    'location','northwest', 'Orientation', 'vertical'); 
legend('boxoff');
title('Abundance may be underestimated at low SNR')


%% signal above noise to noise

[rho p] = corr(FigureS3.abn_data(:,:),FigureS3.eDiffRel(:,:), 'type', 'Pearson', 'rows', 'complete')
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
    end;
end

[rho_all p_all] = corr(reshape(FigureS3.abn_data(:,:),[],1),reshape(FigureS3.eDiffRel(:,:),[],1), 'type', 'Pearson', 'rows', 'complete')

subplot(1,3,2);
for indCond = 1:numel(FigureS3.cond)
    hdot{indCond} = scatter(FigureS3.abn_data(:,indCond),FigureS3.eDiffRel(:,indCond),[],FigureS3.colorm(indCond,:),'filled');
    hold on;
    Fit = polyfit(FigureS3.abn_data(:,indCond),FigureS3.eDiffRel(:,indCond),1);
    plot(FigureS3.abn_data(:,indCond), polyval(Fit,FigureS3.abn_data(:,indCond)), 'Color', FigureS3.colorm(indCond,:));
end; xlim([0 1]);
h_legend3 = legend([hdot{1}, hdot{2}, hdot{3}], ['Load 2: r = ', num2str(rho(1,1),'%.2f'), '; p = ', scientificVals{1,1}],...
    ['Load 4: r = ', num2str(rho(2,2),'%.2f'), '; p = ', scientificVals{2,2}]',...
    ['Load 6: r = ', num2str(rho(3,3),'%.2f'), '; p = ', scientificVals{3,3}]',...
    'location','northwest', 'Orientation', 'vertical');
legend('boxoff');
xlabel('Abundance'); ylabel('Signal-above-noise-to-noise ratio')
title('Effective peak signal is related to abundance')

%% abundance - BGfit

[rho p] = corr(FigureS3.abn_data(:,:),FigureS3.e_BGfit(:,:), 'type', 'Pearson', 'rows', 'complete')
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

[rho_all p_all] = corr(reshape(FigureS3.abn_data(:,:),[],1),reshape(FigureS3.e_BGfit(:,:),[],1), 'type', 'Pearson', 'rows', 'complete')

subplot(1,3,3);
for indCond = 1:numel(FigureS3.cond)
    hdot{indCond} = scatter(FigureS3.abn_data(:,indCond),FigureS3.e_BGfit(:,indCond),[],FigureS3.colorm(indCond,:),'filled');
    hold on;
    Fit = polyfit(FigureS3.abn_data(:,indCond),FigureS3.e_BGfit(:,indCond),1);
    plot(FigureS3.abn_data(:,indCond), polyval(Fit,FigureS3.abn_data(:,indCond)), 'Color', FigureS3.colorm(indCond,:));
end; xlim([0 1]);
h_legend4 = legend([hdot{1}, hdot{2}, hdot{3}], ['Load 2: r = ', num2str(rho(1,1),'%.2f'), '; p = ', scientificVals{1,1}],...
    ['Load 4: r = ', num2str(rho(2,2),'%.2f'), '; p = ', scientificVals{2,2}]',...
    ['Load 6: r = ', num2str(rho(3,3),'%.2f'), '; p = ', scientificVals{3,3}]',...
    'location','northwest', 'Orientation', 'vertical');
legend('boxoff');
xlabel('Abundance'); ylabel('Background amplitude [a.u.]')
title({'Background amplitudes are less';'consistently related to abundance'})

set(findall(gcf,'-property','FontSize'),'FontSize',15)
set(h_legend1, 'FontSize', 12);
set(h_legend3, 'FontSize', 12);
set(h_legend4, 'FontSize', 12);

%% save Figure

pn.plotFolder = ['/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F_Figures/'];
figureName = 'S3';

saveas(h, [pn.plotFolder, figureName], 'fig');
saveas(h, [pn.plotFolder, figureName], 'epsc');
saveas(h, [pn.plotFolder, figureName], 'png');