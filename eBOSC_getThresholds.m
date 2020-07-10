function [eBOSC, pt, dt] = eBOSC_getThresholds(cfg, TFR, e)

    % average power estimates across periods of interest
    BG = [];
    for indTrial = 1:eBOSC.Ntrial
        % remove BGpad at beginning and end to avoid edge artifacts
        BG = [BG TFR.trial{indTrial}(:,cfg.eBOSC.pad.background_sample+1:end-cfg.eBOSC.pad.background_sample)];
    end; clear indTrial
    
    % if frequency ranges should be exluded to reduce the influence of
    % rhythmic peaks on the estimation of the linear background, the
    % following section removes these specified ranges
    freqKeep = true(size(cfg.eBOSC.F));
    if ~isempty(cfg.eBOSC.threshold.excludePeak) % allow for no peak removal
        for indExFreq = 1:size(cfg.eBOSC.threshold.excludePeak,1) % allow for multiple peaks
            % find empirical peak in specified range
            freqInd1 = find(cfg.eBOSC.F >= cfg.eBOSC.threshold.excludePeak(indExFreq,1), 1, 'first');
            freqInd2 = find(cfg.eBOSC.F <= cfg.eBOSC.threshold.excludePeak(indExFreq,2), 1, 'last');
            [~, indPos] = max(mean(BG(freqInd1:freqInd2,:),2));
            indPos = freqInd1+indPos;
            % approximate wavelet extension in frequency domain
            % note: we do not remove the specified range, but the FWHM
            % around the empirical peak
            LowFreq = cfg.eBOSC.F(indPos)-(((2/cfg.eBOSC.wavenumber)*cfg.eBOSC.F(indPos))/2);
            UpFreq = cfg.eBOSC.F(indPos)+(((2/cfg.eBOSC.wavenumber)*cfg.eBOSC.F(indPos))/2);
            % index power estimates within the above range to remove from BG fit
            freqKeep(cfg.eBOSC.F >= LowFreq & cfg.eBOSC.F <= UpFreq) = 0;
        end
    end
    fitInput.f_ = cfg.eBOSC.F(freqKeep);
    fitInput.BG_ = BG(freqKeep, :);
        
    % perform the robust linear fit, only including putatively aperiodic components (i.e., peak exclusion)
    b = robustfit(log10(fitInput.f_),mean(log10(fitInput.BG_),2)'); clear fitInput;
    pv(1) = b(2); pv(2) = b(1);
    mp = 10.^(polyval(pv,log10(cfg.eBOSC.F))); 

    % compute eBOSC power (pt) and duration (dt) thresholds: 
    % power threshold is based on a chi-square distribution with df=2 and mean as estimated above
    pt=chi2inv(cfg.eBOSC.threshold.percentile,2)*mp/2; % chi2inv.m is part of the statistics toolbox of Matlab and Octave
    % duration threshold is the specified number of cycles, so it scales with frequency
    dt=(cfg.eBOSC.threshold.duration*cfg.eBOSC.fsample./cfg.eBOSC.F)';

    % save multiple time-invariant estimates that could be of interest:
    % overall wavelet power spectrum (NOT only background)
    eBOSC.static.bg_pow(e,:)        = mean(BG(:,cfg.eBOSC.pad.total_sample+1:end-cfg.eBOSC.pad.total_sample),2);
    % log10-transformed wavelet power spectrum (NOT only background)
    eBOSC.static.bg_log10_pow(e,:)  = mean(log10(BG(:,cfg.eBOSC.pad.total_sample+1:end-cfg.eBOSC.pad.total_sample)),2);
    % intercept and slope parameters of the robust linear 1/f fit (log-log)
    eBOSC.static.pv(e,:)            = pv;
    % linear background power at each estimated frequency
    eBOSC.static.mp(e,:)            = mp;
    % statistical power threshold
    eBOSC.static.pt(e,:)            = pt;

end
