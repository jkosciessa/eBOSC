%    This file is part of the extended Better OSCillation detection (eBOSC) library.
%
%    The eBOSC library is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    The eBOSC library is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
%
%    Copyright 2020 Julian Q. Kosciessa, Thomas H. Grandy, Douglas D. Garrett & Markus Werkle-Bergner

function [eBOSC, pt, dt] = eBOSC_getThresholds(cfg, TFR, eBOSC)
% This function estimates the static duration and power thresholds and
% saves information regarding the overall spectrum and background.
%
% Inputs: 
%           cfg | config structure with cfg.eBOSC field
%           TFR | time-frequency matrix
%           eBOSC | main eBOSC output structure; will be updated
%
% Outputs: 
%           eBOSC   | updated w.r.t. background info (see below)
%                   | bg_pow: overall power spectrum
%                   | bg_log10_pow: overall power spectrum (log10)
%                   | pv: intercept and slope of fit
%                   | mp: linear background power
%                   | pt: power threshold
%           pt | empirical power threshold
%           dt | duration threshold


    % average power estimates across periods of interest
    BG = [];
    for indTrial = 1:numel(cfg.eBOSC.trial_background)
        % remove BGpad at beginning and end to avoid edge artifacts
        BG = [BG TFR.trial{cfg.eBOSC.trial_background(indTrial)}(:,cfg.eBOSC.pad.background_sample+1:end-cfg.eBOSC.pad.background_sample)];
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
            freqidx = freqInd1:freqInd2;
            [~, indPos] = max(mean(BG(freqidx,:),2));
            indPos = freqidx(indPos);
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
    eBOSC.static.bg_pow(cfg.tmp.channel(1),:)        = mean(BG(:,:),2);
    % log10-transformed wavelet power spectrum (NOT only background)
    eBOSC.static.bg_log10_pow(cfg.tmp.channel(1),:)  = mean(log10(BG(:,:)),2);
    % intercept and slope parameters of the robust linear 1/f fit (log-log)
    eBOSC.static.pv(cfg.tmp.channel(1),:)            = pv;
    % linear background power at each estimated frequency
    eBOSC.static.mp(cfg.tmp.channel(1),:)            = mp;
    % statistical power threshold
    eBOSC.static.pt(cfg.tmp.channel(1),:)            = pt;

end
