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

function [pv,meanpower]=BOSC_bgfit(F,B)
% [pv,meanpower]=BOSC_bgfit(F,B)
%
% This function estimates the background power spectrum via a
% linear regression fit to the power spectrum in log-log coordinates
% 
% parameters:
% F - vector containing frequencies sampled
% B - matrix containing power as a function of frequency (rows) and
% time). This is the time-frequency data.
%
% returns:
% pv = contains the slope and y-intercept of regression line
% meanpower = mean power values at each frequency 
%

pv=polyfit(log10(F),mean(log10(B),2)',1); % linear regression

% transform back to natural units (power; usually uV^2/Hz)
meanpower=10.^(polyval(pv,log10(F)));     
