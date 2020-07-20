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
%    Copyright 2018 Julian Q. Kosciessa, Thomas H. Grandy, Douglas D. Garrett & Markus Werkle-Bergner.

function [pv,meanpower]=eBOSC_bgfit_robust(F,B)
%
% This function estimates the background power spectrum via a
% linear regression fit to the power spectrum in log-log coordinates
% 
% parameters:
% F - vector containing the sampled frequencies (i.e. potentially removing any unwanted frequency peaks)
% B - matrix containing power as a function of the frequencies stated in F (rows) and
% time). This is the time-frequency data.
%
% returns:
% pv = contains the slope and y-intercept of regression line
% meanpower = mean power values at each frequency 
%

b = robustfit(log10(F),mean(log10(B),2)');
pv(1) = b(2);
pv(2) = b(1);

% transform back to natural units (power; usually uV^2/Hz)
meanpower=10.^(polyval(pv,log10(F))); 