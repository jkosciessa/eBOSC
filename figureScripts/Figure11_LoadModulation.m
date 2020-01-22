%% plot topoplots of the Load effect during the last 2 s

% add convertPtoExponential
addpath('/Volumes/EEG/BOSC/BOSC_Sternberg/T_tools/')
addpath('/Volumes/EEG/BOSC/BOSC_Sternberg/T_tools/fieldtrip-20180227'); ft_defaults;

%% plot only theta + alpha

load('/Volumes/EEG/BOSC/BOSC_Sternberg/B_analyses/C2_TFR/B_data/B_CBPA_load_avgAcrossSessions_final2s.mat', 'stat');

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
        title({'Baselined theta (4-8 Hz) power', ['Parametric load; p = ', pval{1}]});
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
        title({'Baselined alpha (8-15 Hz) power', ['Parametric load; p = ', pval{1}]});
        cb = colorbar; set(get(cb,'label'),'string','t values');

%% get abundance etc. results, merge in single figure

parameters = {'o_data', 'a_data', 'e_data', 'o_BGfit', 'a_BGfit', 'e_BGfit', ...
    'abn_data', 'oDiff', 'aDiff', 'eDiff', 'oRel', 'aRel', 'eRel', 'oDiffRel', ...
    'aDiffRel', 'eDiffRel', 'e_iaf', 'e_PT'};

load('/Volumes/EEG/BOSC/BOSC_Sternberg/B_analyses/B_BOSC/B_eBOSC_180917_retention/X_LoadEffect_eBOSCestimates/B_data/X6_stat_avgAcrossSessions_byFreq.mat', 'stat', 'cfgStat')

cfg = [];
cfg.layout = 'acticap-64ch-standard2.mat';
cfg.parameter = 'powspctrm';
cfg.comment = 'no';
cfg.colorbar = 'EastOutside';
cfg.zlim = [-5 5];

% figure;
subplot(2,3,2); cla;
    indParameter = 7; indFreq = 2;
    cfg.highlight = 'yes';
    cfg.highlightchannel = stat{indParameter,indFreq}.label(stat{indParameter,indFreq}.mask(:));
    plotData = [];
    plotData.label = stat{indParameter,indFreq}.label; % {1 x N}
    plotData.dimord = 'chan';
    plotData.powspctrm = squeeze(stat{indParameter,indFreq}.stat(:));
    ft_topoplotER(cfg,plotData);
    cb = colorbar; set(get(cb,'label'),'string','t values');
    pval = []; pval = convertPtoExponential(stat{indParameter,indFreq}.posclusters(1).prob);
    title({'Rhythmic theta abundance', ['Parametric load; p = ', pval{1}]});
subplot(2,3,3+2); cla;
    indParameter = 7; indFreq = 3;
    cfg.highlight = 'yes';
    cfg.highlightchannel = stat{indParameter,indFreq}.label(stat{indParameter,indFreq}.mask(:));
    plotData = [];
    plotData.label = stat{indParameter,indFreq}.label; % {1 x N}
    plotData.dimord = 'chan';
    plotData.powspctrm = squeeze(stat{indParameter,indFreq}.stat(:));
    ft_topoplotER(cfg,plotData);
    cb = colorbar; set(get(cb,'label'),'string','t values');
    pval = []; pval = convertPtoExponential(stat{indParameter}.posclusters(1).prob);
    title({'Rhythmic alpha abundance', ['Parametric load; p = ', pval{1}]});
subplot(2,3,3); cla;
    indParameter = 3; indFreq = 2;
    cfg.highlight = 'yes';
    cfg.highlightchannel = stat{indParameter,indFreq}.label(stat{indParameter,indFreq}.mask(:));
    plotData = [];
    plotData.label = stat{indParameter,indFreq}.label; % {1 x N}
    plotData.dimord = 'chan';
    plotData.powspctrm = squeeze(stat{indParameter,indFreq}.stat(:));
    ft_topoplotER(cfg,plotData);
    cb = colorbar; set(get(cb,'label'),'string','t values');
    title({'Rhythmic theta amplitude'; 'No clusters detected'});
subplot(2,3,3+3); cla;
    indParameter = 3; indFreq = 3;
    cfg.highlight = 'yes';
    cfg.highlightchannel = stat{indParameter,indFreq}.label(stat{indParameter,indFreq}.mask(:));
    plotData = [];
    plotData.label = stat{indParameter,indFreq}.label; % {1 x N}
    plotData.dimord = 'chan';
    plotData.powspctrm = squeeze(stat{indParameter,indFreq}.stat(:));
    ft_topoplotER(cfg,plotData);
    cb = colorbar; set(get(cb,'label'),'string','t values');
    title({'Rhythmic alpha amplitude'; 'No clusters detected'});

% add colorbrewer
addpath('/Volumes/LNDG/Projects/StateSwitch/dynamic/data/eeg/rest/B_analyses/A_MSE_CSD_multiVariant/T_tools/brewermap')
cBrew = brewermap(500,'RdBu');
cBrew = flipud(cBrew);
colormap(cBrew)

%% save Figure
    
set(findall(gcf,'-property','FontSize'),'FontSize',16)
pn.plotFolder = '/Volumes/EEG/BOSC/BOSC_SternRest/X_documentation/B_2018_Manuscript/Figures_resubmit1/F1_Figures/';
figureName = 'F11';
saveas(h, [pn.plotFolder, figureName], 'fig');
saveas(h, [pn.plotFolder, figureName], 'epsc');
saveas(h, [pn.plotFolder, figureName], 'png');