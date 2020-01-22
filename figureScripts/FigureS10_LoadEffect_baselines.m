%% plot topoplots of the Load effect during the last 2 s

% add convertPtoExponential
addpath('/Volumes/EEG/BOSC_Sternberg/T_tools/')
addpath('/Volumes/EEG/BOSC_Sternberg/T_tools/fieldtrip-20180227'); ft_defaults;

%% plot only theta + alpha

load('/Volumes/EEG/BOSC_Sternberg/B_analyses/C2_TFR/B_data/B_CBPA_load_avgAcrossSessions_final2s.mat', 'stat');

parameters = {'powspctrm_load_bl', 'powspctrm_load_bl_ST', ...
    'powspctrm_load_log10', 'powspctrm_load_raw'};
parameterLabels = {'average baseline', 'single trial baseline', ...
    'log10-transformed', 'raw'};

cfg = [];
cfg.layout = 'acticap-64ch-standard2.mat';
cfg.parameter = 'powspctrm';
cfg.comment = 'no';
cfg.colorbar = 'EastOutside';
cfg.zlim = [-5 5];

freqs{1} = stat{1}.freq>4 & stat{1}.freq<8;
freqs{2} = stat{1}.freq>8 & stat{1}.freq<15;
freqs{3} = stat{1}.freq>15 & stat{1}.freq<25;
freqs{4} = stat{1}.freq>25 & stat{1}.freq<64;

freqLabels = {'Theta (4-8 Hz)', 'Alpha (8-15 Hz)'};

h = figure('units','normalized','position',[.1 .1 .5 .3]);
indParameter = 1;
subplot(2,3,1)
        indFreq = 1;
        cfg.highlight = 'yes';
        cfg.highlightchannel = stat{indParameter}.label(nanmax(nanmax(stat{indParameter}.posclusterslabelmat(:,freqs{indFreq},:)==2,[],3),[],2)>0);
        plotData = [];
        plotData.label = stat{indParameter}.label; % {1 x N}
        plotData.dimord = 'chan';
        plotData.powspctrm = squeeze(nanmean(stat{indParameter}.stat(:,freqs{indFreq}),3));
        ft_topoplotER(cfg,plotData);
        title(freqLabels{indFreq})
        cb = colorbar; set(get(cb,'label'),'string','t values');
        pval = []; pval = convertPtoExponential(stat{indParameter}.posclusters(2).prob);
        title({'Baselined Theta (4-8 Hz) power', ['Parametric load; p = ', pval{1}]});
subplot(2,3,3+1)
        indFreq = 2;
        cfg.highlight = 'yes';
        cfg.highlightchannel = stat{indParameter}.label(nanmax(nanmax(stat{indParameter}.posclusterslabelmat(:,freqs{indFreq},:)==3,[],3),[],2)>0);
        plotData = [];
        plotData.label = stat{indParameter}.label; % {1 x N}
        plotData.dimord = 'chan';
        plotData.powspctrm = squeeze(nanmean(stat{indParameter}.stat(:,freqs{indFreq}),3));
        ft_topoplotER(cfg,plotData);
        pval = []; pval = convertPtoExponential(stat{indParameter}.posclusters(3).prob);
        title({'log10 Alpha (8-15 Hz) power', ['Parametric load; p = ', pval{1}]});
        cb = colorbar; set(get(cb,'label'),'string','t values');
indParameter = 3;
subplot(2,3,2)
        indFreq = 1;
        cfg.highlight = 'yes';
        cfg.highlightchannel = stat{indParameter}.label(nanmax(nanmax(stat{indParameter}.posclusterslabelmat(:,freqs{indFreq},:)==2,[],3),[],2)>0);
        plotData = [];
        plotData.label = stat{indParameter}.label; % {1 x N}
        plotData.dimord = 'chan';
        plotData.powspctrm = squeeze(nanmean(stat{indParameter}.stat(:,freqs{indFreq}),3));
        ft_topoplotER(cfg,plotData);
        title(freqLabels{indFreq})
        cb = colorbar; set(get(cb,'label'),'string','t values');
        pval = []; pval = convertPtoExponential(stat{indParameter}.posclusters(2).prob);
        title({'log10 Theta (4-8 Hz) power', ['Parametric load; p = ', pval{1}]});
subplot(2,3,3+2)
        indFreq = 2;
        cfg.highlight = 'yes';
        cfg.highlightchannel = stat{indParameter}.label(nanmax(nanmax(stat{indParameter}.posclusterslabelmat(:,freqs{indFreq},:)==3,[],3),[],2)>0);
        plotData = [];
        plotData.label = stat{indParameter}.label; % {1 x N}
        plotData.dimord = 'chan';
        plotData.powspctrm = squeeze(nanmean(stat{indParameter}.stat(:,freqs{indFreq}),3));
        ft_topoplotER(cfg,plotData);
        pval = []; pval = convertPtoExponential(stat{indParameter}.posclusters(3).prob);
        title({'Baselined Alpha (8-15 Hz) power', ['Parametric load; p = ', pval{1}]});
        cb = colorbar; set(get(cb,'label'),'string','t values');
indParameter = 4;
subplot(2,3,3)
        indFreq = 1;
        cfg.highlight = 'yes';
        cfg.highlightchannel = stat{indParameter}.label(nanmax(nanmax(stat{indParameter}.posclusterslabelmat(:,freqs{indFreq},:)==2,[],3),[],2)>0);
        plotData = [];
        plotData.label = stat{indParameter}.label; % {1 x N}
        plotData.dimord = 'chan';
        plotData.powspctrm = squeeze(nanmean(stat{indParameter}.stat(:,freqs{indFreq}),3));
        ft_topoplotER(cfg,plotData);
        title(freqLabels{indFreq})
        cb = colorbar; set(get(cb,'label'),'string','t values');
        pval = []; pval = convertPtoExponential(stat{indParameter}.posclusters(2).prob);
        title({'Raw Theta (4-8 Hz) power', ['Parametric load; p = ', pval{1}]});
subplot(2,3,3+3)
        indFreq = 2;
        cfg.highlight = 'yes';
        cfg.highlightchannel = stat{indParameter}.label(nanmax(nanmax(stat{indParameter}.posclusterslabelmat(:,freqs{indFreq},:)==4,[],3),[],2)>0);
        plotData = [];
        plotData.label = stat{indParameter}.label; % {1 x N}
        plotData.dimord = 'chan';
        plotData.powspctrm = squeeze(nanmean(stat{indParameter}.stat(:,freqs{indFreq}),3));
        ft_topoplotER(cfg,plotData);
        pval = []; pval = convertPtoExponential(stat{indParameter}.posclusters(4).prob);
        title({'Raw Alpha (8-15 Hz) power', ['Parametric load; p = ', pval{1}]});
        cb = colorbar; set(get(cb,'label'),'string','t values');

% add colorbrewer
addpath('/Volumes/LNDG/Projects/StateSwitch/dynamic/data/eeg/rest/B_analyses/A_MSE_CSD_multiVariant/T_tools/brewermap')
cBrew = brewermap(500,'RdBu');
cBrew = flipud(cBrew);
colormap(cBrew)

%% save Figure
    
set(findall(gcf,'-property','FontSize'),'FontSize',16)
pn.plotFolder = '/Volumes/EEG/BOSC_Sternberg/B_analyses/C2_TFR/C_figures/';
figureName = 'Z_statsOverview_v3_normalization';
saveas(h, [pn.plotFolder, figureName], 'fig');
saveas(h, [pn.plotFolder, figureName], 'epsc');
saveas(h, [pn.plotFolder, figureName], 'png');