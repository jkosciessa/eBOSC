% Reproduce Supplemental Figure 1

%% load measures etc.

% CHOOSE OUTPUT FOLDERS

pn.data = '/Volumes/EEG/BOSC_Simulation/17_JQK_RED/D_data/D_simulation/';
pn.plotFolder = '/Volumes/EEG/BOSC_Simulation/17_JQK_RED/B_figures/';

% Load simulation results

Standard = load([pn.data, 'REDSimulation_standardBOSC_170630.mat']);

ExtendedA = load([pn.data, 'REDSimulation_171023_v10A.mat']);
ExtendedB = load([pn.data, 'REDSimulation_171023_v10B.mat']);
ExtendedC = load([pn.data, 'REDSimulation_171023_v10C.mat']);
ExtendedD = load([pn.data, 'REDSimulation_171023_v10D.mat']);

abundance_ep(1,:,:,:) = ExtendedA.abundance_ep(:,:,:); % rhythmic episode abundance
abundance_ep(2,:,:,:) = ExtendedB.abundance_ep(:,:,:);
abundance_ep(3,:,:,:) = ExtendedC.abundance_ep(:,:,:);
abundance_ep(4,:,:,:) = ExtendedD.abundance_ep(:,:,:);

abundance_PepMeanAlpha = Standard.abundance_PepMeanAlpha; % pepisode of average detected in alpha range
abundance_meanPep = Standard.abundance_meanPep; % average of pepisodes in alpha range

% Load information about simulation

% amount of simulated sample points (according to cycles below)
AmountExt = ExtendedA.Amount;
AmountStd = Standard.Amount;
Amount = AmountExt;

amplitude = [0 2 4 6 8 12 16 24];
cycles = [2 4 8 16 32 64 128 200];
alphaFreq = 10;

amountAlpha = round(round((cycles/alphaFreq),3)/0.004,0);
amountAlpha(end) = 3500; % final abundance is 1, i.e. covering the entire period

% method x amplitude x cycles x repetitions [for extended BOSC]

% correct for wrong indexing on 170302:
% correct size is 1x8x8x200

Standard.SignalDetection.Hits = Standard.SignalDetection.Hits(1,1:8,1:8,:);
Standard.SignalDetection.Misses = Standard.SignalDetection.Misses(1,1:8,1:8,:);
Standard.SignalDetection.CRs = Standard.SignalDetection.CRs(1,1:8,1:8,:);
Standard.SignalDetection.FAs = Standard.SignalDetection.FAs(1,1:8,1:8,:);

abundance_PepMeanAlpha = reshape(abundance_PepMeanAlpha, [1, size(abundance_PepMeanAlpha)]);
abundance_meanPep = reshape(abundance_meanPep, [1, size(abundance_meanPep)]);

% create abundance and cycle labels

for c = 1:numel(cycles)
    tmp_amountAlpha = Amount.Alpha(c);
    tmp_amountNoAlpha = Amount.NoAlpha(c);
    Abundance(1,c) = Amount.Alpha(c)./3500;
    cycLabel{1,c} = [num2str(round(cycles(c),2)), [' ' num2str(round(Abundance(c),2)), '']];
end
cycLabel = cellfun(@(x) strrep(x,' ','\newline'), cycLabel,'UniformOutput',false);

% create labels containing the approximated empirical SNR (overall and episode)
% Note that the SNR will vary depending on the abundance, but this is not reflected here.

SNR = squeeze(Standard.SignalDetection.Amp)./squeeze(Standard.SignalDetection.fitBG);
empiricalSNR = round(max(SNR,[],2),0);
for a = 1:numel(amplitude)
    ampSNR{a} = [num2str(amplitude(a)), ' (', num2str(empiricalSNR(a)),') '];
end

SNRe = squeeze(ExtendedA.SignalDetection.Amp)./squeeze(ExtendedA.SignalDetection.fitBG);
SNRe = squeeze(SNRe(:,:,:));
empiricalSNRe = round(max(SNRe,[],2),0);
for a = 1:numel(amplitude)
    ampSNRe{a} = [num2str(amplitude(a)), ' (', num2str(empiricalSNRe(a)),') '];
end

%% Figure S1A

h = figure('units','normalized','position',[.1 .1 .9 .9]);

tmp_data = squeeze(nanmean(abundance_PepMeanAlpha(1,:,:,:),4));
haxis(1) = subplot(4, 3, [1,4]); bar(tmp_data','edgecolor','none','showbaseline','off')
xlim([.2 8.8])
ylim([0 1.02])
for j = 1:8
    hold on; plot([j-.45 j+.45],[Abundance(j) Abundance(j)],'color','k')
end; clear j
title('Standard BOSC')
%xlabel('Simulated Alpha cycles/abundance')
set(gca,'xtick',[1:8],'xticklabel',{})
ylabel('Detected Alpha abundance')
legend('Amplitude = 0','Amplitude = 2','Amplitude = 4','Amplitude = 6','Amplitude = 8','Amplitude = 12','Amplitude = 16','Amplitude = 24', 'location', 'northwest')
legend('boxoff')

tmp_data = squeeze(nanmean(abundance_ep(2,:,:,:),4));
haxis(2) = subplot(4,3,[2,5]); bar(tmp_data','edgecolor','none','showbaseline','off')
xlim([.2 8.8])
ylim([0 1.02])
for j = 1:8
    hold on; plot([j-.45 j+.45],[Abundance(j) Abundance(j)],'color','k')
end; clear j
title('Extended BOSC')
%xlabel('Simulated Alpha cycles/abundance')
set(gca,'xtick',[1:8],'xticklabel',{})

haxis(3) = subplot(4,3,[3,6]);

% grab positions (1,3 - left,right)
pos_1 = get( haxis(1), 'Position' )
pos_2 = get( haxis(2), 'Position' )
pos_3 = get( haxis(3), 'Position' )

tmp_data = squeeze(nanmean(Standard.SignalDetection.Hits(1,:,:,:),4));
haxis(4) = subplot(4,3,7); bar(tmp_data','edgecolor','none','showbaseline','off')
xlim([.2 8.8])
ylim([0 3700])
for j = 1:8
    hold on; plot([j-.45 j+.45],[amountAlpha(j) amountAlpha(j)],'color','k')
end; clear j
%title('Standard BOSC')
%xlabel('Simulated Alpha cycles/abundance')
set(gca,'xtick',[1:8],'xticklabel',{})
set(gca,'ytick',[0 3500],'yticklabel',{'0'; '3500'});
ylabel('Amount of Hits')
%legend('Amplitude = 0','Amplitude = .25','Amplitude = .50','Amplitude = 1','Amplitude = 2','Amplitude = 4','Amplitude = 8','Amplitude = 16', 'location', 'northwest')
%legend('boxoff')
pos = get( haxis(4), 'Position' )
pos(1) = pos_1(1);
pos(2) = pos(2)+(.001*pos_1(2)) ; % Shift up.
pos(3) = pos_1(3);
pos(4) = pos(4)+(.075*pos_1(2)) ; % Increase height.
set( haxis(4), 'Position', pos) ;

tmp_data = squeeze(nanmean(ExtendedB.SignalDetection.Hits(1,:,:,:),4));
haxis(5) = subplot(4,3,8); bar(tmp_data','edgecolor','none','showbaseline','off')
xlim([.2 8.8])
ylim([0 3700])
for j = 1:8
    hold on; plot([j-.45 j+.45],[amountAlpha(j) amountAlpha(j)],'color','k')
end; clear j
%title('Extended BOSC: FWHM')
%xlabel('Simulated Alpha cycles/abundance')
set(gca,'xtick',[1:8],'xticklabel',{})
set(gca,'ytick',[0 3500],'yticklabel',{'0'; '3500'});
pos = get( haxis(5), 'Position' )
pos(1) = pos_2(1);
pos(2) = pos(2)+(.001*pos_2(2)) ; % Shift up.
pos(3) = pos_2(3);
pos(4) = pos(4)+(.075*pos_2(2)) ; % Increase height.
set( haxis(5), 'Position', pos) ;

haxis(6) = subplot(4,3,9);
pos = get( haxis(6), 'Position' )
pos(1) = pos_3(1);
pos(2) = pos(2)+(.001*pos_3(2)) ; % Shift up.
pos(3) = pos_3(3);
pos(4) = pos(4)+(.075*pos_3(2)) ; % Increase height.
set( haxis(6), 'Position', pos) ;

tmp_data = squeeze(nanmean(Standard.SignalDetection.FAs(1,:,:,:),4));
haxis(7) = subplot(4,3,10); bar(tmp_data','edgecolor','none','showbaseline','off')
xlim([.2 8.8])
ylim([0 3700])
for j = 1:8
    hold on; plot([j-.45 j+.45],[3500-amountAlpha(j) 3500-amountAlpha(j)],'color','k')
end; clear j
%title('Standard BOSC')
xlabel('Alpha cycles/abundance')
set(gca,'xtick',[1:8],'xticklabel',cycLabel)
set(gca,'ytick',[0 3500],'yticklabel',{'0'; '3500'});
ylabel('Amount of FAs')
%legend('Amplitude = 0','Amplitude = .25','Amplitude = .50','Amplitude = 1','Amplitude = 2','Amplitude = 4','Amplitude = 8','Amplitude = 16', 'location', 'northwest')
%legend('boxoff')
pos = get( haxis(7), 'Position' )
pos(1) = pos_1(1);
pos(2) = pos(2)+(2*.001*pos_1(2)) ; % Shift up.
pos(3) = pos_1(3);
pos(4) = pos(4)+(.075*pos_1(2)) ; % Increase height.
set( haxis(7), 'Position', pos) ;

tmp_data = squeeze(nanmean(ExtendedB.SignalDetection.FAs(1,:,:,:),4));
haxis(8) = subplot(4,3,11); bar(tmp_data','edgecolor','none','showbaseline','off')
xlim([.2 8.8])
ylim([0 3700])
for j = 1:8
    hold on; plot([j-.45 j+.45],[3500-amountAlpha(j) 3500-amountAlpha(j)],'color','k')
end; clear j
%title('Extended BOSC: MaxBias')
xlabel('Alpha cycles/abundance')
set(gca,'xtick',[1:8],'xticklabel',cycLabel)
set(gca,'ytick',[0 3500],'yticklabel',{'0'; '3500'});
pos = get( haxis(8), 'Position' )
pos(1) = pos_2(1);
pos(2) = pos(2)+(2*.001*pos_2(2)) ; % Shift up.
pos(3) = pos_2(3);
pos(4) = pos(4)+(.075*pos_2(2)) ; % Increase height.
set( haxis(8), 'Position', pos) ;

haxis(9) = subplot(4,3,12);
pos = get( haxis(9), 'Position' )
pos(1) = pos_3(1);
pos(2) = pos(2)+(2*.001*pos_3(2)) ; % Shift up.
pos(3) = pos_3(3);
pos(4) = pos(4)+(.075*pos_3(2)) ; % Increase height.
set( haxis(9), 'Position', pos) ;
        
set(findall(gcf,'-property','FontSize'),'FontSize',16);

%% Figure S1B

empiricalSNR = round(max(SNR,[],2),0);

h = figure('units','normalized','position',[.1 .1 .6 .9]);
subplot(3,2,1);
    curData = squeeze(nanmean(abundance_PepMeanAlpha(1,:,:,:),4))-repmat((Amount.Alpha./3500),8,1);
    imagesc(curData, [-.2, .2]);
    textStrings = num2str(curData(:),'%0.2f');  %# Create strings from the matrix values
    textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
    [x,y] = meshgrid(1:8, 1:8);   %# Create x and y coordinates for the strings
    hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
                    'HorizontalAlignment','center');
    midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
    textColors = repmat(curData(:) < midValue,1,3);
    textColors(isnan(curData(:)),:) = 1;
    set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
    ylabel('Amplitudes (Empirical SNR)'); xlabel('Cycles (Abundance)');
    set(gca, 'YTick', 1:8);
    set(gca, 'YTickLabels', ampSNR);
    set(gca, 'XTickLabels', cycLabel);
    title('Abundance Error: Standard BOSC');
subplot(3,2,2);
    curData = squeeze(nanmean(abundance_ep(2,:,:,:),4))-repmat((Amount.Alpha./3500),8,1);
    imagesc(curData, [-.2, .2]);
    textStrings = num2str(curData(:),'%0.2f');  %# Create strings from the matrix values
    textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
    [x,y] = meshgrid(1:8, 1:8);   %# Create x and y coordinates for the strings
    hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
                    'HorizontalAlignment','center');
    midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
    textColors = repmat(curData(:) < midValue,1,3);
    textColors(isnan(curData(:)),:) = 1;
    set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
    ylabel('Amplitudes (Empirical SNR)'); xlabel('Cycles (Abundance)');
    set(gca, 'YTick', 1:8);
    set(gca, 'YTickLabels', ampSNR);
    set(gca, 'XTickLabels', cycLabel);
    title('Abundance Error: Extended BOSC');
subplot(3,2,3);
    curData = squeeze(nanmean(Standard.SignalDetection.HitRate(1,:,:,:),4));
    imagesc(curData, [0 1]);
    textStrings = num2str(curData(:),'%0.2f');  %# Create strings from the matrix values
    textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
    [x,y] = meshgrid(1:8, 1:8);   %# Create x and y coordinates for the strings
    hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
                    'HorizontalAlignment','center');
    midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
    textColors = repmat(curData(:) < midValue,1,3);
    textColors(isnan(curData(:)),:) = 1;
    set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
    ylabel('Amplitudes (Empirical SNR)'); xlabel('Cycles (Abundance)');
    set(gca, 'YTick', 1:8);
    set(gca, 'YTickLabels', ampSNR);
    set(gca, 'XTickLabels', cycLabel);
    title('Hit Rate: Standard BOSC');
subplot(3,2,4);
    curData = squeeze(nanmean(ExtendedB.SignalDetection.HitRate(:,:,:),3));
    imagesc(curData, [0 1]);
    textStrings = num2str(curData(:),'%0.2f');  %# Create strings from the matrix values
    textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
    [x,y] = meshgrid(1:8, 1:8);   %# Create x and y coordinates for the strings
    hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
                    'HorizontalAlignment','center');
    midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
    textColors = repmat(curData(:) < midValue,1,3);
    textColors(isnan(curData(:)),:) = 1;
    set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
    ylabel('Amplitudes (Empirical SNR)'); xlabel('Cycles (Abundance)');
    set(gca, 'YTick', 1:8);
    set(gca, 'YTickLabels', ampSNR);
    set(gca, 'XTickLabels', cycLabel);
    title('Hit Rate: ExtendedBOSC');
subplot(3,2,5);
    curData = squeeze(nanmean(Standard.SignalDetection.FARate(1,:,:,:),4));
    imagesc(curData, [0 1]);
    textStrings = num2str(curData(:),'%0.2f');  %# Create strings from the matrix values
    textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
    [x,y] = meshgrid(1:8, 1:8);   %# Create x and y coordinates for the strings
    hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
                    'HorizontalAlignment','center');
    midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
    textColors = repmat(curData(:) < midValue,1,3);
    textColors(isnan(curData(:)),:) = 1;
    set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
    ylabel('Amplitudes (Empirical SNR)'); xlabel('Cycles (Abundance)');
    set(gca, 'YTick', 1:8);
    set(gca, 'YTickLabels', ampSNR);
    set(gca, 'XTickLabels', cycLabel);
    title('FA Rate: Standard BOSC');
subplot(3,2,6);
    curData = squeeze(nanmean(ExtendedB.SignalDetection.FARate(:,:,:),3));
    imagesc(curData, [-.2 .2]);
    textStrings = num2str(curData(:),'%0.2f');  %# Create strings from the matrix values
    textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
    [x,y] = meshgrid(1:8, 1:8);   %# Create x and y coordinates for the strings
    hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
                    'HorizontalAlignment','center');
    midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
    textColors = repmat(curData(:) < midValue,1,3);
    textColors(isnan(curData(:)),:) = 1;
    set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
    ylabel('Amplitudes (Empirical SNR)'); xlabel('Cycles (Abundance)');
    set(gca, 'YTick', 1:8);
    set(gca, 'YTickLabels', ampSNR);
    set(gca, 'XTickLabels', cycLabel);
    title('FA Rate: Extended BOSC');
    
set(findall(gcf,'-property','FontSize'),'FontSize',15)
% saveas(h, [pn.plotFolder, datestr(now, 'yymmdd'), '_RelativePlot_HM'], 'fig');
% saveas(h, [pn.plotFolder, datestr(now, 'yymmdd'), '_RelativePlot_HM'], 'epsc');
