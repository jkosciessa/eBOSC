% Plot eBOSC detection example

% - requires B_BOSC/A_detectionExample/A_scripts/A_CalculateExampleBOSCDetection_180605.m
% - original at /Volumes/EEG/BOSC_Sternberg/scripts/B_BOSC/A_detectionExample/A_scripts/B_plotDetectionExample.m

% load('/Volumes/EEG/BOSC_Sternberg/B_analyses/B_BOSC/A_detectionExample/B_data/DetectionExample.mat', 'dataToPlot')
% Figure1BCDE = dataToPlot;
% 
% addpath('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F4_FigureTools/brewermap/')
% Figure1BCDE.colorm = brewermap(4,'RdYlBu'); % use external function to create colormap
% 
% %% save Figure data
% 
% save('/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F3_FigureData/F1BCDE.mat', 'Figure1BCDE');

%% load Figure data

load('/Users/kosciessa/Desktop/eBOSC/figureData/F1BCDE.mat', 'Figure1BCDE');

%% plot 

h = figure('units','normalized','position',[.1 .1 .5 .7]);
subplot(2,2,1)
plot(Figure1BCDE.param.BOSC.F, 10.^mean(log10(Figure1BCDE.BG),2), 'k', 'LineWidth', 2); hold on;
plot(Figure1BCDE.param.BOSC.F, 10.^polyval(Figure1BCDE.pv_std,log10(Figure1BCDE.param.BOSC.F)), 'Color', Figure1BCDE.colorm(4,:), 'LineWidth', 2)
plot(Figure1BCDE.param.BOSC.F, 10.^polyval(Figure1BCDE.pv,log10(Figure1BCDE.param.BOSC.F)), 'Color', Figure1BCDE.colorm(1,:), 'LineWidth', 2)
xlim([0 64]);
title('Background fit'); xlabel('Frequency [Hz]'); ylabel('Power [a.u.]');
subplot(2,2,2)
plot(log10(Figure1BCDE.param.BOSC.F), mean(log10(Figure1BCDE.BG),2), 'k', 'LineWidth', 2); hold on;
plot(log10(Figure1BCDE.param.BOSC.F), polyval(Figure1BCDE.pv_std,log10(Figure1BCDE.param.BOSC.F)), 'Color', Figure1BCDE.colorm(4,:), 'LineWidth', 2)
plot(log10(Figure1BCDE.param.BOSC.F), polyval(Figure1BCDE.pv,log10(Figure1BCDE.param.BOSC.F)), 'Color', Figure1BCDE.colorm(1,:), 'LineWidth', 2)
title('Background fit (log-log)'); xlabel('log10 Frequency [log10 Hz]'); ylabel('log10 Power [a.u.]');
legend({'original signal'; 'linear fit'; 'robust fit'}, 'location', 'southwest');
legend('boxoff');

tmp.axisTicks = [5 10 15 20 25 30 35 40];
h1 = subplot(2,2,3);
ax = gca;
imagesc(Figure1BCDE.B_(:,Figure1BCDE.param.BOSC.trialPad+1:end-Figure1BCDE.param.BOSC.trialPad)); hold on;
originalSize = get(gca, 'Position');
visboundaries(ax, Figure1BCDE.detected(:,Figure1BCDE.param.BOSC.detectedPad+1:end-Figure1BCDE.param.BOSC.detectedPad), 'LineWidth', 0.5, 'Color', 'w');
visboundaries(ax, Figure1BCDE.detected_abn(:,Figure1BCDE.param.BOSC.detectedPad+1:end-Figure1BCDE.param.BOSC.detectedPad), 'LineWidth', 3, 'Color', 'r');
title('Power with detected segments');
set(gca,'YDir','normal'); 
hcol = colorbar;
ylabel(hcol, 'Power [a.u.]')
ax = gca;
ax.YTick = tmp.axisTicks;
ax.YTickLabel = round(Figure1BCDE.param.BOSC.F(tmp.axisTicks),1);
ax.XTick = [150 375 600];
ax.XTickLabel = [Figure1BCDE.time(ax.XTick)*1000];
xlabel('Time [ms]'); ylabel('Frequency [Hz]');
caxis([0*10^5 5*10^5])
v = caxis;

% calculate the abundance within the alpha range

tmp.idx = cell2mat(Figure1BCDE.episodes(:,7))>= 8 & cell2mat(Figure1BCDE.episodes(:,7))<= 13; % Original segments
matchEpSPs = cell2mat(Figure1BCDE.episodes(tmp.idx,5)); % Original segments
uniqueSPs = unique(matchEpSPs(:,2));
abundance_alpha = numel(uniqueSPs)/750;

alpha_idx = find(Figure1BCDE.param.BOSC.F >= 8 & Figure1BCDE.param.BOSC.F <= 13);

subplot(2,2,4);
tempplot = repmat(Figure1BCDE.pt',1,size(Figure1BCDE.B_(:,Figure1BCDE.param.BOSC.trialPad+1:end-Figure1BCDE.param.BOSC.trialPad),2));
for indRow = 1:numel(Figure1BCDE.dt)
    tempplot(indRow,ceil(Figure1BCDE.dt(indRow)):end) = NaN;
end;
imagesc(tempplot, v);
ax = gca;
visboundaries(ax, ~isnan(tempplot), 'LineWidth', 3, 'Color', 'k');
set(gca,'YDir','normal');
hold on; scatter(Figure1BCDE.pepisode*size(tempplot,2), 1:numel(Figure1BCDE.pt), 400,'w.');
line([abundance_alpha*size(tempplot,2) abundance_alpha*size(tempplot,2)], [min(alpha_idx) max(alpha_idx)], 'Color', [1,0,0], 'LineWidth', 10)
ax = gca; ax.YTick = []; 
ax.XTick = [150 375 600];
ax.XTickLabel = round([ax.XTick/750],3);
caxis([0*10^5 5*10^5])

title('Power and Duration threshold');
xlabel('Abundance'); ylabel('Frequency (see left)')
axesHandles = findobj(get(h,'Children'), 'flat','Type','axes');
axis(axesHandles,'square')
set(h1, 'Position', originalSize);

set(findall(gcf,'-property','FontSize'),'FontSize',19)

colormap(parula)

pn.plotFolder = '/Volumes/EEG/BOSC_SternRest/X_documentation/B_2018_Manuscript/F_Figures/';
figureName = 'F1BCDE';

saveas(h, [pn.plotFolder, figureName], 'fig');
saveas(h, [pn.plotFolder, figureName], 'epsc');
saveas(h, [pn.plotFolder, figureName], 'png');