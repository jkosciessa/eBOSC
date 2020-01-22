%% Plot Figure F9D (triadic split of rhythmic amplitudes)

pn.dataOut = '/Volumes/EEG/BOSC_Sternberg/B_analyses/B_BOSC/B_eBOSC_180509_full_wl6_noDur/B_180509_MergeRhythmCharacteristics/B_data/';
load([pn.dataOut, 'D_eventsSplitByeBOSC.mat'], 'smallEvents_wavelet*', 'mediumEvents_wavelet*', 'largeEvents_wavelet*')

time = -500:500; time = time.*.004;

%% Plot rhythm representation (eAmp) with statistics as inlay

% add within-subject error bars
pn.shadedError = ['/Volumes/EEG/BOSC_Sternberg/T_tools/shadedErrorBar-7002ebc']; addpath(pn.shadedError);

h = figure('units','normalized','position',[0 0 .3 .4]);
cla; hold on;

events_eAmp = cat(4, smallEvents_wavelet_e_amp, mediumEvents_wavelet_e_amp,largeEvents_wavelet_e_amp);

    condAvg = squeeze(nanmean(nanmean(events_eAmp,2),4));
    curData = squeeze(nanmean(nanmean(events_eAmp(:,:,:,1),2),4));
    curData = curData-condAvg+repmat(nanmean(condAvg,1),32,1);
    standError = nanstd(curData,1)./sqrt(size(curData,1));
    l1 = shadedErrorBar(time*1000,nanmean(curData,1),standError, 'lineprops', {'color', [.6 .6 .6],'linewidth', 3}, 'patchSaturation', .05);
    condAvg = squeeze(nanmean(nanmean(events_eAmp,2),4));
    curData = squeeze(nanmean(nanmean(events_eAmp(:,:,:,2),2),4));
    curData = curData-condAvg+repmat(nanmean(condAvg,1),32,1);
    standError = nanstd(curData,1)./sqrt(size(curData,1));
    l2 = shadedErrorBar(time*1000,nanmean(curData,1),standError, 'lineprops', {'color', [.3 .3 .3],'linewidth', 3}, 'patchSaturation', .05);
    condAvg = squeeze(nanmean(nanmean(events_eAmp,2),4));
    curData = squeeze(nanmean(nanmean(events_eAmp(:,:,:,3),2),4));
    curData = curData-condAvg+repmat(nanmean(condAvg,1),32,1);
    standError = nanstd(curData,1)./sqrt(size(curData,1));
    l3 = shadedErrorBar(time*1000,nanmean(curData,1),standError, 'lineprops', {'color', [0 0 0],'linewidth', 3}, 'patchSaturation', .05);
xlabel('Time in ms'); xlim([-500, 500]); ylabel('Amplitude (?Volt)');ylim([-17 10])
title({'Single-trial rhythmic amplitude estimates';'reflect variations in time series amplitude'})
ll1 = legend([l1.mainLine, l2.mainLine, l3.mainLine], ...
    {'Small amplitude'; 'Medium amplitude'; 'Large amplitude'}, 'location', 'SouthWest'); legend('boxoff');

% add stats inlay

addpath(['/Volumes/Kosciessa/Tools/barwitherr/']);
addpath(['/Volumes/Kosciessa/Tools/mysigstar/']);

idxTime = find(time >= -.1 & time < .1);

events_eAmp = cat(4, smallEvents_wavelet_e_amp, mediumEvents_wavelet_e_amp,largeEvents_wavelet_e_amp).^2;
events_oAmp = cat(4, smallEvents_wavelet_o_amp, mediumEvents_wavelet_o_amp,largeEvents_wavelet_o_amp).^2;
events_eAbn = cat(4, smallEvents_wavelet_e_abn, mediumEvents_wavelet_e_abn,largeEvents_wavelet_e_abn).^2;

condPairs = [1,2; 2,3; 1,3];
condPairsLevel = [50 65 80];

handaxes2 = axes('Position', [0.68 0.25 0.2 .2]);

    curData = squeeze(nanmean(nanmean(events_eAmp(:,1:3,idxTime,:),3),2))';% avg over time, conditions
    meanY = nanmean(curData,2);
    condAvg = squeeze(nanmean(curData,1)); % avg over conditions
    curData_within = curData-repmat(condAvg,3,1)+repmat(nanmean(condAvg,2),3,32);
        
    errorY = nanstd(curData_within,[],2)/sqrt(size(curData_within,2));
    [h1, hError] = barwitherr(errorY, meanY);
    for indPair = 1:size(condPairs,1)
        % significance star for the difference
        [~, pval] = ttest(curData(condPairs(indPair,1),:), curData(condPairs(indPair,2), :)); % paired t-test
        % if mysigstar gets 2 xpos inputs, it will draw a line between them and the
        % sigstars on top
        PvalPairs(indPair) = pval;
        if pval <.05
            mysigstar(gca, [condPairs(indPair,1)+.1 condPairs(indPair,2)-.1], condPairsLevel(indPair), pval);
        end
    end
    set(h1(1),'FaceColor',[.8 .3 .25]);
    set(h1(1),'LineWidth',2);
    set(hError(1),'LineWidth',3);
    box(gca,'off')
    ylim([0 80])
    set(gca, 'XTickLabels', {'S', 'M', 'L'})
    xlabel({'Amplitude bin'}); ylabel({'Power';'[?Volt^2]'})

    set(findall(gcf,'-property','FontSize'),'FontSize',23)

pn.plotFolder = '/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/Figures_resubmit1/F1_Figures/';
figureName = 'Figure9D_triadicSplit';

saveas(h, [pn.plotFolder, figureName], 'fig');
saveas(h, [pn.plotFolder, figureName], 'epsc');
saveas(h, [pn.plotFolder, figureName], 'png');