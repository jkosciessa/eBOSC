%    This file is part of the Better OSCillation detection (BOSC) library.
%
%    The BOSC library is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    The BOSC library is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
%
%    Copyright 2010 Jeremy B. Caplan, Adam M. Hughes, Tara A. Whitten
%    and Clayton T. Dickson.

function [powthresh,durthresh]=BOSC_thresholds(Fsample,percentilethresh,numcyclesthresh,F,meanpower)
% 
% This function calculates all the power thresholds and duration
% thresholds for use with BOSC_detect.m to detect oscillatory episodes
% Fsample - sampling rate (Hz)
% percentilethresh - power threshold expressed as a percentile/100
%                    (i.e., from 0-1) of the estimated
%                    chi-square(2) probability distribution of
%                    power values. A typical value is 0.95
% numcyclesthresh - duration threshold. A typical value is 3 cycles.
% 
% returns:
% power thresholds and duration thresholds
%
% Fsample = is the sampling rate
%
% F - frequencies sampled in the power spectrum
%
% meanpower - power spectrum (mean power at each frequency)


% power threshold is based on a chi-square distribution with df=2
% and mean as estimated previously (BOSC_bgfit.m)
powthresh=chi2inv(percentilethresh,2)*meanpower/2;
% chi2inv.m is part of the statistics toolbox of Matlab and Octave

% duration threshold is simply a certain number of cycles, so it
% scales with frequency
durthresh=(numcyclesthresh*Fsample./F)';
