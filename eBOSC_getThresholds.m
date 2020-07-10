function [eBOSC] = eBOSC_getThresholds(cfg, TFR)

    % average power estimates across periods of interest
    BG = [];
    for indTrial = 1:eBOSC.Ntrial
        % remove BGpad at beginning and end to avoid edge artifacts
        BG = [BG TFR.trial{indTrial}(:,cfg.eBOSC.pad.background_sample+1:end-cfg.eBOSC.pad.background_sample)];
    end; clear indTrial

    % To DO: allow for NO peak removal, multiple peak removal

    cfg.eBOSC.static.excludePeak = [5, 8; 20 30];

    cfg.eBOSC.LowFreqExcludeBG = 

    % find peak in specified range
    freqInd1 = find(cfg.eBOSC.F >= cfg.eBOSC.LowFreqExcludeBG, 1, 'first');
    freqInd2 = find(cfg.eBOSC.F <= cfg.eBOSC.HighFreqExcludeBG, 1, 'last');
    [~, indPos] = max(mean(BG(freqInd1:freqInd2,:),2));
    indPos = freqInd1+indPos;
    % approximate wavelet extension in frequency domain
    LowFreq = cfg.eBOSC.F(indPos)-(((2/cfg.eBOSC.wavenumber)*cfg.eBOSC.F(indPos))/2);
    UpFreq = cfg.eBOSC.F(indPos)+(((2/cfg.eBOSC.wavenumber)*cfg.eBOSC.F(indPos))/2);
    % remove power estimates within the following range from the fit
    freqIndLow = find(cfg.eBOSC.F >= LowFreq, 1, 'first');
    freqIndHigh = find(cfg.eBOSC.F <= UpFreq, 1, 'last');
    % perform the robust linear fit, only including putatively aperiodic components (i.e., peak exclusion)
    [pv,~] = eBOSC_bgfit_robust(cfg.eBOSC.F([1:freqIndLow-1 freqIndHigh+1:end]),BG([1:freqIndLow-1 freqIndHigh+1:end], :));
    mp = 10.^(polyval(pv,log10(cfg.eBOSC.F))); 

    % compute eBOSC power (pt) and duration (dt) thresholds: 
    [pt,dt] = BOSC_thresholds(cfg.eBOSC.fsample,cfg.eBOSC.percentile,cfg.eBOSC.ncyc,cfg.eBOSC.F,mp);

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
