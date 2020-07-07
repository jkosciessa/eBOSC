load('/Volumes/EEG/BOSC/BOSC_Sternberg/B_analyses/B_BOSC/B_eBOSC_180509_full_wl6_noDur/B_180509_MergeRhythmCharacteristics/B_data/D_stat.mat', 'stat_burst', 'stat_rhythm', 'cfgStat')

%% plot load effects

% add colorbrewer
addpath('/Volumes/LNDG/Projects/StateSwitch/dynamic/data/eeg/rest/B_analyses/A_MSE_CSD_multiVariant/T_tools/brewermap')
cBrew = brewermap(500,'RdBu');
cBrew = flipud(cBrew);

cfg = [];
cfg.layout = 'acticap-64ch-standard2.mat';
cfg.parameter = 'powspctrm';
cfg.comment = 'no';
cfg.colorbar = 'EastOutside';
cfg.style = 'both';
cfg.colormap = cBrew;
cfg.zlim = [-5 5];


paramNames = {'# of episodes', 'Cycle duration', 'Frequency', 'Power'};
selectParams = [7,1,3,5]; locations = [4,8,12,16];
h = figure('units','normalized','position',[0 0 1 1]);
for indParameter = 1:numel(selectParams)
    curParam = selectParams(indParameter);
    indFreq = 1;
    subplot(4,4,locations(indParameter));
    cfg.highlight = 'yes';
    cfg.highlightchannel = stat_burst{curParam,indFreq}.label(stat_burst{curParam,indFreq}.mask);
    plotData = [];
    plotData.label = stat_burst{curParam,indFreq}.label; % {1 x N}
    plotData.dimord = 'chan';
    plotData.powspctrm = squeeze(stat_burst{curParam,indFreq}.stat);
    if numel(find(isnan(plotData.powspctrm)))==60
        continue;
    else
        ft_topoplotER(cfg,plotData);
    end
    cb = colorbar; set(get(cb,'label'),'string','t values');
    if indParameter == 1
        pval = []; pval = convertPtoExponential(stat_burst{curParam,indFreq}.posclusters(1).prob);
        title({[paramNames{indParameter}], ['probe vs. retention: p = ', pval{1}]});
    elseif indParameter == 2
        pval = []; pval = convertPtoExponential(stat_burst{curParam,indFreq}.negclusters(1).prob);
        title({[paramNames{indParameter}], ['probe vs. retention: p = ', pval{1}]});
    elseif indParameter == 3
        pval = []; pval = convertPtoExponential(stat_burst{curParam,indFreq}.negclusters(1).prob);
        title({[paramNames{indParameter}], ['probe vs. retention: p = ', pval{1}]});
    elseif indParameter == 4
        pval = []; pval = convertPtoExponential(stat_burst{curParam,indFreq}.negclusters(1).prob);
        title({[paramNames{indParameter}], ['probe vs. retention: p = ', pval{1}]});
    end
end

paramNames = {'# of episodes', 'Cycle duration', 'Frequency', 'Power'};
selectParams = [7,1,3,5]; locations = [1,5,9,13];
%h = figure('units','normalized','position',[0 0 .25 1]);
for indParameter = 1:numel(selectParams)
    curParam = selectParams(indParameter);
    indFreq = 1;
    subplot(4,4,locations(indParameter));
    cfg.highlight = 'yes';
    cfg.highlightchannel = stat_rhythm{curParam,indFreq}.label(stat_rhythm{curParam,indFreq}.mask);
    plotData = [];
    plotData.label = stat_rhythm{curParam,indFreq}.label; % {1 x N}
    plotData.dimord = 'chan';
    plotData.powspctrm = squeeze(stat_rhythm{curParam,indFreq}.stat);
    if numel(find(isnan(plotData.powspctrm)))==60
        continue;
    else
        ft_topoplotER(cfg,plotData);
    end
    cb = colorbar; set(get(cb,'label'),'string','t values');
    if indParameter == 1
        pval = []; pval = convertPtoExponential(stat_rhythm{curParam,indFreq}.posclusters(1).prob);
        title({[paramNames{indParameter}], ['retention vs. probe: p = ', pval{1}]});
    elseif indParameter == 2
        pval = []; pval = convertPtoExponential(stat_rhythm{curParam,indFreq}.posclusters(1).prob);
        title({[paramNames{indParameter}], ['retention vs. probe: p = ', pval{1}]});
    elseif indParameter == 4
        pval = []; pval = convertPtoExponential(stat_rhythm{curParam,indFreq}.posclusters(1).prob);
        title({[paramNames{indParameter}], ['retention vs. probe: p = ', pval{1}]});
    else
        title([paramNames{indParameter}])
    end
end

%% load Figure data

addpath('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F4_FigureTools/');

load('/Users/kosciessa/Desktop/eBOSC/figureData/F10.mat', 'Figure14')
addpath('/Volumes/EEG/BOSC/BOSC_Sternberg/T_tools') % needs convertPtoExponential

%% Figure: Rhythm vs. burst characteristics

%h = figure('units','normalized','position',[.1 .1 .7 .9]); 

colormap(cBrew)
subplot(4,4,2);imagesc(Figure14.timevec, [], Figure14.DetSumRhythm, [0 50]); title('Median Detected Episodes');
cb1 = colorbar; set(get(cb1,'label'),'string','# of episodes');
addTaskTiming(get(gca,'YLim')); xlabel('Time (s)'); ylabel({'Channel';'(anterior-posterior)'}); xlim([-2 5])
subplot(4,4,6);imagesc(Figure14.timevec, [], Figure14.CycleDurationMeanRhythm); title('Median Cycle Duration');
cb2 = colorbar; set(get(cb2,'label'),'string','Duration (cycles)');
addTaskTiming(get(gca,'YLim')); xlabel('Time (s)'); ylabel({'Channel';'(anterior-posterior)'}); xlim([-2 5])
subplot(4,4,10);imagesc(Figure14.timevec, [], Figure14.FreqMeanRhythm, [10 12]); title('Median Frequency');
cb3 = colorbar; set(get(cb3,'label'),'string','Frequency (Hz)');
addTaskTiming(get(gca,'YLim')); xlabel('Time (s)'); ylabel({'Channel';'(anterior-posterior)'}); xlim([-2 5])
subplot(4,4,14);imagesc(Figure14.timevec, [], Figure14.PowMeanRhythm); title('Median Power');
cb4 = colorbar; set(get(cb4,'label'),'string','Power (a.u.)');
addTaskTiming(get(gca,'YLim')); xlabel('Time (s)'); ylabel({'Channel';'(anterior-posterior)'}); xlim([-2 5])

subplot(4,4,3);imagesc(Figure14.timevec, [], Figure14.DetSumBurst, [0 50]); title('Median Detected Episodes');
cb5 = colorbar; set(get(cb5,'label'),'string','# of episodes');
addTaskTiming(get(gca,'YLim')); xlabel('Time (s)'); ylabel({'Channel';'(anterior-posterior)'}); xlim([-2 5])
subplot(4,4,7);imagesc(Figure14.timevec, [], Figure14.CycleDurationMeanBurst); title('Median Cycle Duration');
cb6 = colorbar; set(get(cb6,'label'),'string','Duration (cycles)');
addTaskTiming(get(gca,'YLim')); xlabel('Time (s)'); ylabel({'Channel';'(anterior-posterior)'}); xlim([-2 5])
subplot(4,4,11);imagesc(Figure14.timevec, [], Figure14.FreqMeanBurst, [10 12]); title('Median Frequency');
cb7 = colorbar; set(get(cb7,'label'),'string','Frequency (Hz)');
addTaskTiming(get(gca,'YLim')); xlabel('Time (s)'); ylabel({'Channel';'(anterior-posterior)'}); xlim([-2 5])
subplot(4,4,15);imagesc(Figure14.timevec, [], Figure14.PowMeanBurst); title('Median Power');
cb8 = colorbar; set(get(cb8,'label'),'string','Power (a.u.)');
addTaskTiming(get(gca,'YLim')); xlabel('Time (s)'); ylabel({'Channel';'(anterior-posterior)'}); xlim([-2 5])

set(findall(gcf,'-property','FontSize'),'FontSize',18)

pn.plotFolder = ['/Volumes/EEG/BOSC/BOSC_SternRest/X_documentation/B_2018_Manuscript/F_Figures/'];
figureName = 'F14_v3_A';

saveas(h, [pn.plotFolder, figureName], 'fig');
saveas(h, [pn.plotFolder, figureName], 'epsc');
saveas(h, [pn.plotFolder, figureName], 'png');

%% add time-series representation

load('/Users/kosciessa/Desktop/eBOSC/figureData/F10C.mat', 'Figure14C')

addpath('/Volumes/EEG/BOSC/BOSC_Sternberg/T_tools/tight_subplot')

h = figure('units','normalized','position',[0 0 1 .25]);

    ax1 = subplot(1,5,[1,2]); cla
        colormap(ax1,'gray');
        imagesc(Figure14C.timeConcat_rhythm(2:end)*1000,[],Figure14C.smoothedVersion_rhythm, [-1 1])
        xlim([-400, 400]); ylabel('Transient episodes')
        ax1.YAxis.Exponent = 4;
        yyaxis right
        hold on;
        for indID = 1:32
            plot(Figure14C.time(2:end)*1000,squeeze(nanmean(Figure14C.AvgSeriesRhythms(:,indID,:),1)), 'Color', [1 .9 0], 'LineWidth', 1, 'Marker', 'none', 'LineStyle', '-')
            xlabel('Time in ms'); xlim([-400, 400]); ylabel('Amplitude (z-scored)')
        end
        hold on; plot(Figure14C.time(2:end)*1000,squeeze(nanmean(nanmean(Figure14C.AvgSeriesRhythms,2),1)), 'Color', [1 0 0], 'LineWidth', 3, 'Marker', 'none', 'LineStyle', '-')
        title('Alpha rhythms during retention (8-15 Hz; O2)')
        ylim([-2.5 2.5])
    ax2 = subplot(1,5,[4,5]); cla;
        colormap(ax2,'gray');
        imagesc(Figure14C.timeConcat_burst(2:end)*1000,[],Figure14C.smoothedVersion_burst, [-1 1])
        xlim([-400, 400]); ylabel('Rhythmic episodes')
        ax2.YAxis.Exponent = 4;
        yyaxis right
        hold on;
        for indID = 1:32
            plot(Figure14C.time(2:end)*1000,squeeze(nanmean(Figure14C.AvgSeriesBurst(:,indID,:),1)), 'Color', [1 .9 0], 'LineWidth', 1, 'Marker', 'none', 'LineStyle', '-')
            xlabel('Time in ms'); xlim([-400, 400]); ylabel('Amplitude (z-scored)')
        end
        hold on; plot(Figure14C.time(2:end)*1000,squeeze(nanmean(nanmean(Figure14C.AvgSeriesBurst,2),1)), 'Color', [1 0 0], 'LineWidth', 3, 'Marker', 'none', 'LineStyle', '-')
        title('Alpha transients during probe (8-15 Hz; O2)')
        ylim([-2.5 2.5])

 subplot(1,5,3); cla; hold on;
    timeIdx_Probe = Figure14.timevec>3 & Figure14.timevec<3.2;
    timeIdx_Retention = Figure14.timevec>0 & Figure14.timevec<3;
    x = squeeze(nanmedian(Figure14.eDiffRel,2)); y = squeeze(nanmedian(Figure14.IndividualDetSum_R(:,timeIdx_Retention),2));
    scatter(x, y, 80,'filled','MarkerFaceColor', Figure14.colorm(1,:));  y1_ls = polyval(polyfit(x,y,1),x); y1_ls = plot(x, y1_ls, 'Color', Figure14.colorm(1,:), 'LineWidth', 2);
    [rl1,pl1] = corrcoef(x,y); pval1 = []; pval1 = convertPtoExponential(pl1(2));
    x = squeeze(nanmedian(Figure14.eDiffRel,2)); y = squeeze(nanmedian(Figure14.IndividualDetSum_A(:,timeIdx_Probe),2));
    scatter(x, y, 80,'filled','MarkerFaceColor', Figure14.colorm(4,:));  y2_ls = polyval(polyfit(x,y,1),x); y2_ls = plot(x, y2_ls, 'Color', Figure14.colorm(4,:), 'LineWidth', 2);
    [rl2,pl2] = corrcoef(x,y); pval2 = []; pval2 = convertPtoExponential(pl2(2));
    leg1 = legend([y1_ls, y2_ls], {['# Detected rhythms:\newliner=',num2str(round(rl1(2),2)), ', p=' pval1{1}];...
            ['# Detected transients:\newliner=',num2str(round(rl2(2),2)), ', p=' pval2{1}]}, 'location', 'NorthWest'); legend('boxoff');
    xlabel('Rhythmic SNR'); ylabel('Amount of detected events')
    ylim([-30 300])
    title({'SNR relates to amount of rhythmic,';' but not transient events'})
    set(findall(gcf,'-property','FontSize'),'FontSize',23)
    set(leg1, 'FontSize', 23)

pn.plotFolder = ['/Volumes/EEG/BOSC/BOSC_SternRest/X_documentation/B_2018_Manuscript/F_Figures/'];
figureName = 'F14C';

saveas(h, [pn.plotFolder, figureName], 'fig');
saveas(h, [pn.plotFolder, figureName], 'epsc');
saveas(h, [pn.plotFolder, figureName], 'png');
   
% set(findall(gcf,'-property','FontSize'),'FontSize',18)
        
%%% Figure: plot interindividual differences in cylce duration sustained vs. burst

% h = figure('units','normalized','position',[.1 .1 .7 .3]); colormap(hot)
% subplot(1,3,1);imagesc(Figure14.timevec, [], Figure14.DetSumRhythm_sorted); title('Median count of rhythmic episodes');
% addTaskTiming(get(gca,'YLim')); xlabel('Time (s)'); ylabel('Subjects (sorted by SNR)'); colorbar; xlim([-2 5])
% 
% subplot(1,3,2);imagesc(Figure14.timevec, [], Figure14.DetSumBurst_sorted,[0 70]); title('Median count of transients episodes');
% addTaskTiming(get(gca,'YLim')); xlabel('Time (s)'); ylabel('Subjects (sorted by SNR)'); colorbar; xlim([-2 5])
% set(findall(gcf,'-property','FontSize'),'FontSize',18)
% 
% subplot(1,3,3); hold on;
% x = squeeze(nanmedian(Figure14.eDiffRel,2)); y = squeeze(nanmedian(Figure14.IndividualDetSum_R,2));
% scatter(x, y, 'filled','MarkerFaceColor', Figure14.colorm(1,:));  y1_ls = polyval(polyfit(x,y,1),x); y1_ls = plot(x, y1_ls, 'Color', Figure14.colorm(1,:), 'LineWidth', 2);
% [rl1,pl1] = corrcoef(x,y); pval1 = []; pval1 = convertPtoExponential(pl1(2));
% x = squeeze(nanmedian(Figure14.eDiffRel,2)); y = squeeze(nanmedian(Figure14.IndividualDetSum_A,2));
% scatter(x, y, 'filled','MarkerFaceColor', Figure14.colorm(4,:));  y2_ls = polyval(polyfit(x,y,1),x); y2_ls = plot(x, y2_ls, 'Color', Figure14.colorm(4,:), 'LineWidth', 2);
% [rl2,pl2] = corrcoef(x,y); pval2 = []; pval2 = convertPtoExponential(pl2(2));
% leg1 = legend([y1_ls, y2_ls], {['# Detected rhythms: r=',num2str(round(rl1(2),2)), ', p=' pval1{1}];...
%         ['# Detected transients: r=',num2str(round(rl2(2),2)), ', p=' pval2{1}]}, 'location', 'NorthWest'); legend('boxoff');
% xlabel('Rhythmic SNR'); ylabel('Amount of detected events')
% title('Rhythm, but not transient rate relates to rhythmic SNR')
% set(findall(gcf,'-property','FontSize'),'FontSize',15)
% set(leg1, 'FontSize', 13)
