% Simulate time series with superimposed rhythm, schematically show
% detected rhythmicity and amplitude gain.

%     pn.root = '/Volumes/EEG/';
%     
% %%  "set path" for analyses
% 
%     % EEG tools folder
%     pn.fun = [pn.root 'ConMemEEGTools/'];
% 
%     % change directory
%     cd([pn.fun 'fnct_THG'])
% 
%     % add function folders folders for analysis
%     THG_addpath_BOSC(pn.fun,0)
%     THG_addpath_FAlpha(pn.fun)
%     THG_addpath_fnct_common(pn.fun)
%     THG_addpath_fnct_thg(pn.fun)
% 
% %%  generate alpha
% 
%     time   = [0:.001:1];
%     alpha1 = sin(2*pi*time*10)*2;
%     
% %%  add 1/f noise
% 
%     % background
%     randn('seed',150901);
%     Figure1A.signal = f_alpha_gaussian(18000,1,1);
%     
%     % add alpha
%     Figure1A.signal(8500:9500) = Figure1A.signal(8500:9500) + alpha1';
%     
%     % filter Figure1A.signal
%     Figure1A.signal = cm_filt_but_lp(Figure1A.signal,1000,100,4);
%     Figure1A.signal = cm_filt_but_hp(Figure1A.signal,1000,.5,4);
% 
% %%  standard BOSC
%     
%     bosc.F           = 2.^[1:.0625:5.25];
%     bosc.fsample     = 1000;
%     bosc.wavenumber  = 6;
%     bosc.padding     = 4;    
%     bosc.percentile  = .99;
%     bosc.numcycles   = 3;
% 
%     [Pepisode,detected,Figure1A.B] = THG_BOSC_standard_wrapper_20150603(Figure1A.signal,bosc);
% 
% %%  new abundance detection
% 
%     abn2.F       = bosc.F;
%     abn2.fsample = bosc.fsample;
%     abn2.fres    = 1.3;
%     abn2.fstp    = 1;
%     abn2.ncyc    = bosc.numcycles;
% 
%     [Figure1A.detected_eBOSC,episodes] = THG_abundance_detect_20150827(Figure1A.B,detected,abn2);
% 
% %% save Figure data
% 
%     save('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/F1A.mat', 'Figure1A');

%% load Figure data

load('/Users/kosciessa/Desktop/eBOSC/figureData/F1A.mat', 'Figure1A');

%%  plot

Figure1A.timevec = [0:0.001:2];

h = figure('units','normalized','position',[.1 .1 .45 .7]);
subplot(4,4,1:3)
    plot(Figure1A.timevec,Figure1A.signal(4000+4000:4000+6000), 'k', 'LineWidth', 2)
    ylim([-6 6])
    ylabel('Amplitude [\muV]')
    box off
subplot(4,4,5:7)
    imagesc(Figure1A.timevec,[], Figure1A.B(:,4000:6000))
    set(gca,'ydir','Normal','ytick',[20 40 60],'yticklabel',{'5','11','26'});%,'xtick',[4000.5 5000 6000],'xticklabel',{'0','1','2'})
    ylabel('Frequency [Hz]')
subplot(4,4,8)
    plot(mean(Figure1A.B(:,4000:6000),2), 1:size(Figure1A.B,1), 'k', 'LineWidth', 2);
    ylim([1 size(Figure1A.B,1)]); xlim([0 5*10^5]);
    set(gca,'ydir','Normal','ytick',[],'yticklabel',{}, 'xtick',{});
    set(gca,'Visible','off')
subplot(4,4,9:11)
    imagesc(Figure1A.timevec,[], Figure1A.B(:,4000:6000))
    [x y] = find(Figure1A.detected_eBOSC(:,4000:6000)==1);
    for j = 1:length(x)
        hold on; plot(Figure1A.timevec(y(j)),x(j),'.r', 'LineWidth', 5)
    end
    set(gca,'ydir','Normal','ytick',[20 40 60],'yticklabel',{'5','11','26'});
    ylabel('Frequency [Hz]')
subplot(4,4,12)
    tmp_B = Figure1A.B;
    tmp_B(Figure1A.detected_eBOSC==0) = NaN;
    tmp_B_range = nanmean(tmp_B(:,4000:6000),2); clear tmp_B;
    tmp_B_range(isnan(tmp_B_range)) = 0;
    plot(tmp_B_range, 1:size(Figure1A.B,1), 'k', 'LineWidth', 2);
    ylim([1 size(Figure1A.B,1)]); xlim([0 5*10^5]);
    set(gca,'ydir','Normal','ytick',[],'yticklabel',{}, 'xtick',{}); xlabel('Power (a.u.)')
    set(gca,'Visible','off')
subplot(4,4,13:15)
    plot(Figure1A.timevec,Figure1A.signal(8000:10000), 'k', 'LineWidth', 2);
    hold on;plot(Figure1A.timevec(400:1629),Figure1A.signal(8399:9628),'r', 'LineWidth', 2);
    ylim([-6 6])
    xlabel('Time [sec]')
    ylabel('Amplitude [\muV]')
    box off

    set(findall(gcf,'-property','FontSize'),'FontSize',18)
    
    colormap(parula)
    
pn.plotFolder = '/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F_Figures/';
figureName = 'F1A';

saveas(h, [pn.plotFolder, figureName], 'fig');
saveas(h, [pn.plotFolder, figureName], 'epsc');
saveas(h, [pn.plotFolder, figureName], 'png');
