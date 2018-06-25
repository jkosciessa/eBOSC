h = figure('units','normalized','position',[.1 .1 .7 .6]);

%% amplitude-abundance correlations

for indVariant = 1:3

    corrAcrossCycles = [];
    for a = 1:numel(amplitude)
        signalOfInterest_Amp = squeeze(data{indVariant}.SignalDetection.BGdiff);
        signalOfInterest_Amp = squeeze(signalOfInterest_Amp(a,:,:));
        signalOfInterest_Abn = squeeze(data{indVariant}.SignalDetection.abn);
        signalOfInterest_Abn = squeeze(signalOfInterest_Abn(a,:,:));
        output = corr(signalOfInterest_Amp(:), signalOfInterest_Abn(:),'rows', 'pairwise');
        corrAcrossCycles(a) = output;
    end

    corrAcrossAmps = [];
    for c = 1:numel(cycles)
        signalOfInterest_Amp = squeeze(data{indVariant}.SignalDetection.BGdiff);
        signalOfInterest_Amp = squeeze(signalOfInterest_Amp(:,c,:));
        signalOfInterest_Abn = squeeze(data{indVariant}.SignalDetection.abn);
        signalOfInterest_Abn = squeeze(signalOfInterest_Abn(:,c,:));
        output = corr(signalOfInterest_Amp(:), signalOfInterest_Abn(:),'rows', 'pairwise');
        corrAcrossAmps(c) = output;
    end

    subplot(2,2,1); 
    hold on; h1(indVariant) = plot(corrAcrossCycles, 'LineWidth', 2); xlabel('Cycles'); ylabel('Correlation coefficient'); ylim([-1 1]);
    hold on; line([0 8], [0 0], 'LineStyle', '--', 'Color', [0 0 0], 'LineWidth', 2)
    set(gca,'xtick',1:8,'xticklabel',cycLabel);
    subplot(2,2,2); 
    hold on; h2(indVariant) = plot(corrAcrossAmps, 'LineWidth', 2); xlabel('Amplitude'); ylabel('Correlation coefficient'); ylim([-1 1]);
    hold on; line([0 8], [0 0], 'LineStyle', '--', 'Color', [0 0 0], 'LineWidth', 2)
    set(gca,'xtick',1:8,'xticklabel',num2cell(amplitude));
end
subplot(2,2,1); 
legend([h1(1), h1(2), h1(3)],{'Standard BOSC'; 'MaxBiasAll'; 'MaxBiasPT'}, 'location', 'south'); legend('boxoff');
title({'Rhythmic Amplitude - abundance';'correlations across amplitudes'});
subplot(2,2,2); 
legend([h2(1), h2(2), h2(3)],{'Standard BOSC'; 'MaxBiasAll'; 'MaxBiasPT'}, 'location', 'south'); legend('boxoff');
title({'Rhythmic Amplitude - abundance';'correlations across cycles'});
set(findall(gcf,'-property','FontSize'),'FontSize',15)

%% abundance-background correlations

for indVariant = 1:3

    corrAcrossCycles = [];
    for a = 1:numel(amplitude)
        signalOfInterest_Amp = squeeze(data{indVariant}.SignalDetection.fitBG);
        signalOfInterest_Amp = squeeze(signalOfInterest_Amp(a,:,:));
        signalOfInterest_Abn = squeeze(data{indVariant}.SignalDetection.abn);
        signalOfInterest_Abn = squeeze(signalOfInterest_Abn(a,:,:));
        output = corr(signalOfInterest_Amp(:), signalOfInterest_Abn(:),'rows', 'pairwise');
        corrAcrossCycles(a) = output;
    end

    corrAcrossAmps = [];
    for c = 1:numel(cycles)
        signalOfInterest_Amp = squeeze(data{indVariant}.SignalDetection.fitBG);
        signalOfInterest_Amp = squeeze(signalOfInterest_Amp(:,c,:));
        signalOfInterest_Abn = squeeze(data{indVariant}.SignalDetection.abn);
        signalOfInterest_Abn = squeeze(signalOfInterest_Abn(:,c,:));
        output = corr(signalOfInterest_Amp(:), signalOfInterest_Abn(:),'rows', 'pairwise');
        corrAcrossAmps(c) = output;
    end

    subplot(2,2,3); 
    hold on; h1(indVariant) = plot(corrAcrossCycles, 'LineWidth', 2); xlabel('Cycles'); ylabel('Correlation coefficient'); ylim([-1 1]);
    hold on; line([0 8], [0 0], 'LineStyle', '--', 'Color', [0 0 0], 'LineWidth', 2)
    set(gca,'xtick',1:8,'xticklabel',cycLabel);
    subplot(2,2,4); 
    hold on; h2(indVariant) = plot(corrAcrossAmps, 'LineWidth', 2); xlabel('Amplitude'); ylabel('Correlation coefficient'); ylim([-1 1]);
    hold on; line([0 8], [0 0], 'LineStyle', '--', 'Color', [0 0 0], 'LineWidth', 2)
    set(gca,'xtick',1:8,'xticklabel',num2cell(amplitude));
end
subplot(2,2,3); 
legend([h1(1), h1(2), h1(3)],{'Standard BOSC'; 'MaxBiasAll'; 'MaxBiasPT'}, 'location', 'south'); legend('boxoff');
title({'Background - abundance';'correlations across amplitudes'});
subplot(2,2,4); 
legend([h2(1), h2(2), h2(3)],{'Standard BOSC'; 'MaxBiasAll'; 'MaxBiasPT'}, 'location', 'south'); legend('boxoff');
title({'Background - abundance';'correlations across cycles'});
set(findall(gcf,'-property','FontSize'),'FontSize',15)

% pn.plotFolder = '/Volumes/EEG/BOSC_Simulation/17_JQK_RED/B_figures/';
% figureName = 'SimCorrelationBGrelAbn_170929_summary';

pn.plotFolder = '/Volumes/EEG/BOSC_SternRest/F_Figures/';
figureName = 'S4';

saveas(h, [pn.plotFolder, figureName], 'fig');
saveas(h, [pn.plotFolder, figureName], 'epsc');
saveas(h, [pn.plotFolder, figureName], 'png');
