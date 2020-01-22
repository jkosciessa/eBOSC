     
clear all; clc;

%% plot bandpassed versions in overview figure: theta, alpha, beta, gamma

  h = figure('units','normalized','position',[.1 .1 .8 .25]);
  ax = subplot(1,2,2); cla
        colormap(ax,'gray');
        load('/Volumes/EEG/BOSC/BOSC_Sternberg/B_analyses/B_BOSC/B_eBOSC_180509_full_wl6_noDur/B_180509_MergeRhythmCharacteristics/B_data/C_TimeDomainRhythmicity_AlphaBurst_lowpass_concat.mat','TimeEps', 'tempSeries*');
        [~, sortIdx] = sort(tempSeriesDurPost, 'ascend');
        tempSeries(isnan(tempSeries)) = 0;
        smoothedVersion = reshape(smooth(tempSeries(sortIdx,:),150, 'moving'),size(tempSeries));
        timeConcat = -1500*1/250:1/250:(-1500*1/250)+size(smoothedVersion,2)*1/250;
        Figure14C.smoothedVersion_burst = smoothedVersion;
        Figure14C.timeConcat_burst = timeConcat;
        imagesc(Figure14C.timeConcat_burst(2:end)*1000,[],Figure14C.smoothedVersion_burst, [-1 1])
        xlim([-400, 400]); ylabel('Rhythmic episodes')
        ax.YAxis.Exponent = 4;
        yyaxis right
        hold on;
        load('/Volumes/EEG/BOSC/BOSC_Sternberg/B_analyses/B_BOSC/B_eBOSC_180509_full_wl6_noDur/B_180509_MergeRhythmCharacteristics/B_data/C_TimeDomainRhythmicity_AlphaBurst_lowpass.mat','TimeEps');
        for indLoad = 1:3        
        for indID = 1:32
               curData = cat(1, TimeEps{indID,indLoad});
               ampData = curData(:,4);
               curDataPlot = curData(:,3);
               curData = curData(:,1);
               tempSeries = NaN(size(curData,1),3000);
               for indEp = 1:size(curData,1)
                     % find minima clostest to TFR maximum 
                    [~, sortIndTFR] = max(ampData{indEp});
                    TF = islocalmin(curData{indEp});
                    TF = find(TF);
                    [~, minInd_tmp] = min(abs(sortIndTFR-TF));
                    sortInd = TF(minInd_tmp);
                    tempSeries(indEp, 1500-sortInd+1:1500+numel(curData{indEp})-sortInd) = curDataPlot{indEp};
                    tempSeriesDur(indEp) = numel(curData{indEp});
                    tempSeriesDurPre(indEp) = 500+sortInd;
                    tempSeriesDurPost(indEp) = numel(curData{indEp})-500-sortInd;
               end
               AvgSeries(indLoad,indID,:) = nanmedian(tempSeries,1);
               time = -1500*1/250:1/250:1500*1/250;
        end
        end
        Figure14C.AvgSeriesBurst = AvgSeries;
        Figure14C.time = time;
        for indID = 1:32
            plot(Figure14C.time(2:end)*1000,squeeze(nanmean(Figure14C.AvgSeriesBurst(:,indID,:),1)), 'Color', [1 .9 0], 'LineWidth', 1, 'Marker', 'none', 'LineStyle', '-')
            xlabel('Time in ms'); xlim([-400, 400]); ylabel('Amplitude (z-scored)')
        end
        hold on; plot(Figure14C.time(2:end)*1000,squeeze(nanmean(nanmean(Figure14C.AvgSeriesBurst,2),1)), 'Color', [1 0 0], 'LineWidth', 3, 'Marker', 'none', 'LineStyle', '-')
        title('Alpha transients during probe (8-15 Hz; O2)')
        ylim([-2.5 2.5])
    
    ax = subplot(1,2,1); cla
        colormap(ax,'gray');
        tempSeriesDurPost = [];
        load('/Volumes/EEG/BOSC/BOSC_Sternberg/B_analyses/B_BOSC/B_eBOSC_180509_full_wl6_noDur/B_180509_MergeRhythmCharacteristics/B_data/C_TimeDomainRhythmicity_AlphaRhythm_Retention_lowpass_concat.mat','TimeEps', 'tempSeries*');
        [~, sortIdx] = sort(tempSeriesDurPost, 'ascend');
        tempSeries(isnan(tempSeries)) = 0;
        smoothedVersion = reshape(smooth(tempSeries(sortIdx,:),150, 'moving'),size(tempSeries));
        timeConcat = -1500*1/250:1/250:(-1500*1/250)+size(smoothedVersion,2)*1/250;
        Figure14C.smoothedVersion_rhythm = smoothedVersion;
        Figure14C.timeConcat_rhythm = timeConcat;
        imagesc(Figure14C.timeConcat_rhythm(2:end)*1000,[],Figure14C.smoothedVersion_rhythm, [-1 1])
        imagesc(timeConcat(2:end)*1000,[],smoothedVersion, [-1 1])
        xlim([-400, 400]); ylabel('Transient episodes')
        ax.YAxis.Exponent = 4;
        yyaxis right
        hold on;
        load('/Volumes/EEG/BOSC/BOSC_Sternberg/B_analyses/B_BOSC/B_eBOSC_180509_full_wl6_noDur/B_180509_MergeRhythmCharacteristics/B_data/C_TimeDomainRhythmicity_AlphaRyhthm_retention_lowpass.mat','TimeEps');
        for indLoad = 1:3        
        for indID = 1:32
               curData = cat(1, TimeEps{indID,indLoad});
               ampData = curData(:,4);
               curDataPlot = curData(:,3);
               curData = curData(:,1);
               tempSeries = NaN(size(curData,1),3000);
               for indEp = 1:size(curData,1)
                     % find minima clostest to TFR maximum 
                    [~, sortIndTFR] = max(ampData{indEp});
                    TF = islocalmin(curData{indEp});
                    TF = find(TF);
                    [~, minInd_tmp] = min(abs(sortIndTFR-TF));
                    sortInd = TF(minInd_tmp);
                    tempSeries(indEp, 1500-sortInd+1:1500+numel(curData{indEp})-sortInd) = curDataPlot{indEp};
                    tempSeriesDur(indEp) = numel(curData{indEp});
                    tempSeriesDurPre(indEp) = 500+sortInd;
                    tempSeriesDurPost(indEp) = numel(curData{indEp})-500-sortInd;
               end
               AvgSeries(indLoad,indID,:) = nanmedian(tempSeries,1);
               time = -1500*1/250:1/250:1500*1/250;
        end
        end
        Figure14C.AvgSeriesRhythms = AvgSeries;
        Figure14C.time = time;
        for indID = 1:32
            plot(Figure14C.time(2:end)*1000,squeeze(nanmean(Figure14C.AvgSeriesRhythms(:,indID,:),1)), 'Color', [1 .9 0], 'LineWidth', 1, 'Marker', 'none', 'LineStyle', '-')
            xlabel('Time in ms'); xlim([-400, 400]); ylabel('Amplitude (z-scored)')
        end
        hold on; plot(Figure14C.time(2:end)*1000,squeeze(nanmean(nanmean(Figure14C.AvgSeriesRhythms,2),1)), 'Color', [1 0 0], 'LineWidth', 3, 'Marker', 'none', 'LineStyle', '-')
        title('Alpha rhythms during retention (8-15 Hz; O2)')
        ylim([-2.5 2.5])

       set(findall(gcf,'-property','FontSize'),'FontSize',18)

    pn.plotFolder = '/Volumes/EEG/BOSC/BOSC_Sternberg/B_analyses/B_BOSC/B_eBOSC_180509_full_wl6_noDur/B_180509_MergeRhythmCharacteristics/C_figures/';
    figureName = 'C_timeSeriesRhythms_ProbeTransients';

    saveas(h, [pn.plotFolder, figureName], 'fig');
    saveas(h, [pn.plotFolder, figureName], 'epsc');
    saveas(h, [pn.plotFolder, figureName], 'png');

    
    save('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/F14C.mat', 'Figure14C')