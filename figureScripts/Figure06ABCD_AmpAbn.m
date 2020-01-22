% - used to be /Volumes/FB-LIP/ConMem/EEG/BOSC_SternRest/B_scripts_JQK/B_processing_170704/E_JQKplotAmpAbn_170816B_Fig8.m

% %% get data
% 
% setup.dataVersion = '170703';
% setup.scriptDir = '170704';
% setup.plotVersion = '170816B';
% setup.savePlot = 0;
% 
% addpath('/Volumes/EEG/BOSC_SternRest/B_analyses/A_eBOSC/A_SternBRest_170703/C_interrelations/T_tools/brewermap/')
% Figure7.colorm = brewermap(4,'RdYlBu'); % use external function to create colormap
% 
% % load simulated data
% 
% % load(['/Volumes/EEG/BOSC_Simulation/17_JQK_RED/D_data/D_simulation/_SimDataforRestingStateComparison_v4.mat'], 'simData');
% % Figure7.simData = simData;
% Figure7.overallSim = load(['/Volumes/EEG/BOSC_Simulation/17_JQK_RED/D_data/D_simulation/_SimDataforRestingStateComparison_overall_v4.mat'], 'simData');
% 
% info.channels = 44:60;
% 
% addpath(genpath('/Volumes/EEG/BOSC_Sternberg/scripts/Tools/'))
% 
% load(['/Volumes/EEG/BOSC_SternRest/B_analyses/A_eBOSC/A_SternBRest_170703/B_extractIndices/B_data/X15B_8to15_',setup.dataVersion,'.mat']);
% 
% cond = {'EC1','EC2','EO1','EO2'};
% 
% for c = 1:numel(cond)
%     Figure7.e_BGfit(:,c) = nanmean(nanmean(X{1,1}.([cond{c},'e_fitBG'])(:,info.channels,:),3),2)';
%     Figure7.abn_data(:,c) = nanmean(nanmean(X{1,1}.([cond{c},'e_abn'])(:,info.channels,:),3),2)';
%     % lagged coherence
%     Figure7.lC_AMax_data(:,c) = nanmean(nanmean(X{1,1}.([cond{c},'o_lC_AMax'])(:,info.channels,:),3),2)';
%     % xAmp./BG
%     Figure7.oRel(:,c) = nanmean(nanmean(X{1,1}.([cond{c},'o_amp'])(:,info.channels,:)./X{1,1}.([cond{c},'o_fitBG'])(:,info.channels,:),3),2)';
%     % xAmp-BG./BG
%     Figure7.eDiffRel(:,c) = nanmean(nanmean((X{1,1}.([cond{c},'e_amp'])(:,info.channels,:)-X{1,1}.([cond{c},'e_fitBG'])(:,info.channels,:))./X{1,1}.([cond{c},'e_fitBG'])(:,info.channels,:),3),2)';
%     % Pepisode
%     Figure7.Pep(:,c) = nanmean(nanmean(X{1,1}.([cond{c},'o_Pep_orig_IAF'])(:,info.channels,:),3),2)'; % Pepisode @ IAF
% end;

%% save Figure data

%save('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/F7.mat', 'Figure7')

%% load Figure data

% load('/Users/kosciessa/Desktop/eBOSC/figureData/F6.mat', 'Figure7')

%% Plot

h = figure('units','normalized','position',[.1 .1 .7 .7],'DefaultAxesFontSize',15);

%% SNR (overall, NOT rhythm-exclusive)

% Calculate separate correlations fits below an empirical SNR of 6 and above (for all data here).

% below 6

tmp_abn = reshape(Figure7.abn_data, [],1);
tmp_snr = reshape(Figure7.oRel, [],1);

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

subplot(2,2,1)
for indCond = 1:4
    scatter(Figure7.abn_data(:,indCond),(Figure7.oRel(:,indCond)),[],Figure7.colorm(indCond,:),'filled');
    hold on;
end
hold on; plot(Figure7.overallSim.simData.Abn, Figure7.overallSim.simData.SNR, 'k', 'lineWidth', 1.5);
Fit = polyfit(tmp_abn_bel,tmp_snr_bel,1);
h1 = plot(tmp_abn_bel, polyval(Fit,tmp_abn_bel));
Fit = polyfit(tmp_abn_abv,tmp_snr_abv,1);
h2 = plot(tmp_abn_abv, polyval(Fit,tmp_abn_abv));
xlabel('Abundance'); ylabel('Overall SNR')
h_legend1 = legend([h2, h1], ['SNR above 6: r = ', num2str(rho_abv,'%.2f'), '; p = ', scientificVals_abv{1,1}], ...
    ['SNR below 6: r = ', num2str(rho_bel,'%.2f'), '; p = ', scientificVals_bel{1,1}]', ...
    'location','northwest', 'Orientation', 'vertical'); 
legend('boxoff');
title('Abundance may be underestimated at low SNR')


%% lagged coherence

[rho p] = corr(Figure7.abn_data(:,:),(Figure7.lC_AMax_data(:,:)), 'type', 'Pearson', 'rows', 'complete')

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
% ! maximum lagged coherence in alpha range

subplot(2,2,4);
for indCond = 1:4
    hdot{indCond} = scatter(Figure7.abn_data(:,indCond),Figure7.lC_AMax_data(:,indCond),[],Figure7.colorm(indCond,:),'filled');
    hold on;
    Fit = polyfit(Figure7.abn_data(:,indCond),Figure7.lC_AMax_data(:,indCond),1);
    plot(Figure7.abn_data(:,indCond), polyval(Fit,Figure7.abn_data(:,indCond)), 'Color', Figure7.colorm(indCond,:));
end; ylim([0 1]); xlim([0 1]);
xlabel('Abundance'); ylabel('Lagged Coherence (max. 8-15 Hz)')
set(findall(gcf,'-property','FontSize'),'FontSize',18)
h_legend2 = legend([hdot{1}, hdot{2}, hdot{3}, hdot{4}], ['EC: r = ', num2str(rho(1,1),'%.2f'), '; p = ', scientificVals{1,1}],...
    ['ECEO - EC: r = ', num2str(rho(2,2),'%.2f'), '; p = ', scientificVals{2,2}]',...
    ['EO: r = ', num2str(rho(3,3),'%.2f'), '; p = ', scientificVals{3,3}]',...
    ['ECEO - EO: r = ', num2str(rho(4,4),'%.2f'), '; p = ', scientificVals{4,4}]',...
    'location','northwest', 'Orientation', 'vertical');
legend('boxoff');
title({'Abundance captures rhythmicity';'similarly to phase-based estimate'})

%% rhythmic SNR

[rho p] = corr(Figure7.abn_data(:,:),Figure7.eDiffRel(:,:), 'type', 'Pearson', 'rows', 'complete')
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

[rho_all p_all] = corr(reshape(Figure7.abn_data(:,:),[],1),reshape(Figure7.eDiffRel(:,:),[],1), 'type', 'Pearson', 'rows', 'complete')

subplot(2,2,2);
for indCond = 1:4
    hdot{indCond} = scatter(Figure7.abn_data(:,indCond),Figure7.eDiffRel(:,indCond),[],Figure7.colorm(indCond,:),'filled');
    hold on;
    Fit = polyfit(Figure7.abn_data(:,indCond),Figure7.eDiffRel(:,indCond),1);
    plot(Figure7.abn_data(:,indCond), polyval(Fit,Figure7.abn_data(:,indCond)), 'Color', Figure7.colorm(indCond,:));
end; xlim([0 1]);
h_legend3 = legend([hdot{1}, hdot{2}, hdot{3}, hdot{4}], ['EC: r = ', num2str(rho(1,1),'%.2f'), '; p = ', scientificVals{1,1}],...
    ['ECEO - EC: r = ', num2str(rho(2,2),'%.2f'), '; p = ', scientificVals{2,2}]',...
    ['EO: r = ', num2str(rho(3,3),'%.2f'), '; p = ', scientificVals{3,3}]',...
    ['ECEO - EO: r = ', num2str(rho(4,4),'%.2f'), '; p = ', scientificVals{4,4}]',...
    'location','northwest', 'Orientation', 'vertical');
legend('boxoff');
xlabel('Abundance'); ylabel('Rhythmic SNR')
title('Effective peak signal is related to abundance')

%% abundance - BGfit

[rho p] = corr(Figure7.abn_data(:,:),Figure7.e_BGfit(:,:), 'type', 'Pearson', 'rows', 'complete')
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

[rho_all p_all] = corr(reshape(Figure7.abn_data(:,:),[],1),reshape(Figure7.e_BGfit(:,:),[],1), 'type', 'Pearson', 'rows', 'complete')

subplot(2,2,3);
for indCond = 1:4
    hdot{indCond} = scatter(Figure7.abn_data(:,indCond),Figure7.e_BGfit(:,indCond),[],Figure7.colorm(indCond,:),'filled');
    hold on;
    Fit = polyfit(Figure7.abn_data(:,indCond),Figure7.e_BGfit(:,indCond),1);
    plot(Figure7.abn_data(:,indCond), polyval(Fit,Figure7.abn_data(:,indCond)), 'Color', Figure7.colorm(indCond,:));
end; xlim([0 1]);
h_legend4 = legend([hdot{1}, hdot{2}, hdot{3}, hdot{4}], ['EC: r = ', num2str(rho(1,1),'%.2f'), '; p = ', scientificVals{1,1}],...
    ['ECEO - EC: r = ', num2str(rho(2,2),'%.2f'), '; p = ', scientificVals{2,2}]',...
    ['EO: r = ', num2str(rho(3,3),'%.2f'), '; p = ', scientificVals{3,3}]',...
    ['ECEO - EO: r = ', num2str(rho(4,4),'%.2f'), '; p = ', scientificVals{4,4}]',...
    'location','northwest', 'Orientation', 'vertical');
legend('boxoff');
xlabel('Abundance'); ylabel('Background amplitude [a.u.]')
title({'Background amplitudes are not';'consistently related to abundance'})

set(findall(gcf,'-property','FontSize'),'FontSize',17)
set(h_legend1, 'FontSize', 12);
set(h_legend2, 'FontSize', 12);
set(h_legend3, 'FontSize', 12);
set(h_legend4, 'FontSize', 12);

%% save plot

pn.plotFolder = '/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F_Figures/';

saveas(h, [pn.plotFolder, 'F7'], 'fig');
saveas(h, [pn.plotFolder, 'F7'], 'epsc');
saveas(h, [pn.plotFolder, 'F7'], 'png');

%% calculate Pepisode-abundance correlation to report in Manuscript
% Pepisode is based on the trial-wise IAF peak

[rho p] = corr(Figure7.abn_data(:,:),(Figure7.Pep(:,:)), 'type', 'Pearson', 'rows', 'complete')
