% original: /Volumes/EEG/BOSC_Simulation/17_JQK_RED/A_scripts/G_SupplementFigure_v2.m
% similar to /Volumes/EEG/BOSC_Simulation/17_JQK_RED/A_scripts/G_ResultsPlots_extVSstnd_Paper_v10_171023.m

%% load measures etc.

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

%AmountExt = load([pn.data, 'Amount_170626.mat'], 'Amount'); % vector with amount of simulated alpha sample points
AmountExt = ExtendedA.Amount;

%AmountStd = load([pn.data, 'Amount_standardBOSC_170410.mat'], 'Amount'); % vector with amount of simulated alpha sample points
AmountStd = Standard.Amount;

% The above should be the same.

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

%% create overview plot

empiricalSNR = round(max(SNR,[],2),0);

h = figure('units','normalized','position',[.1 .1 .9 .9]);
subplot(5,3,2);
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
subplot(5,3,5);
    curData = squeeze(nanmean(ExtendedA.SignalDetection.HitRate(:,:,:),3));
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
    title('Hit Rate: MaxBias raw');
subplot(5,3,8);
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
    title('Hit Rate: MaxBias PT');
subplot(5,3,3);
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
subplot(5,3,6);
    curData = squeeze(nanmean(ExtendedA.SignalDetection.FARate(:,:,:),3));
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
    title('FA Rate: MaxBias raw');
subplot(5,3,9);
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
    title('FA Rate: MaxBias PT');
subplot(5,3,1);
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
subplot(5,3,4);
    curData = squeeze(nanmean(abundance_ep(1,:,:,:),4))-repmat((Amount.Alpha./3500),8,1);
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
    title('Abundance Error: MaxBias raw');
subplot(5,3,7);
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
    title('Abundance Error: MaxBias PT');
subplot(5,3,11);
    curData = squeeze(nanmean(ExtendedC.SignalDetection.HitRate(:,:,:),3));
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
    title('Hit Rate: FWHM raw');
subplot(5,3,14);
    curData = squeeze(nanmean(ExtendedD.SignalDetection.HitRate(:,:,:),3));
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
    title('Hit Rate: FWHM PT');
subplot(5,3,12);
    curData = squeeze(nanmean(ExtendedC.SignalDetection.FARate(:,:,:),3));
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
    title('FA Rate: FWHM raw');
subplot(5,3,15);
    curData = squeeze(nanmean(ExtendedD.SignalDetection.FARate(:,:,:),3));
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
    title('FA Rate: FWHM PT');
subplot(5,3,10);
    curData = squeeze(nanmean(abundance_ep(3,:,:,:),4))-repmat((Amount.Alpha./3500),8,1);
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
    title('Abundance Error: FWHM raw');
subplot(5,3,13);
    curData = squeeze(nanmean(abundance_ep(4,:,:,:),4))-repmat((Amount.Alpha./3500),8,1);
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
    title('Abundance Error: FWHM PT');