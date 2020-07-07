function [map,scheme] = brewermap_view(N,scheme)
% An interactive figure for ColorBrewer colormap selection. With demo!
%
% (c) 2015 Stephen Cobeldick
%
% View Cynthia Brewer's ColorBrewer colorschemes in a figure.
%
% * Two colorbars give the colorscheme in color and grayscale.
% * A button toggles between 3D-cube and 2D-lineplot of the RGB values.
% * A button toggles an endless cycle through the colorschemes.
% * A button reverses the colormap.
% * 35 buttons select any ColorBrewer colorscheme.
% * Text with the colorscheme's type (Diverging/Qualitative/Sequential)
% * Text with the colorscheme's number of nodes (defining colors).
%
% Syntax:
%  brewermap_view
%  brewermap_view(N)
%  brewermap_view(N,scheme)
%  brewermap_view([],...)
%  brewermap_view({axes/figure handles},...) % see "Adjust External Colormaps"
%  [map,scheme] = brewermap_view(...)
%
% Calling the function with an output argument blocks MATLAB execution until
% the figure is deleted: the final colormap and parameters are then returned.
%
% See also BREWERMAP CUBEHELIX RGBPLOT COLORMAP COLORMAPEDITOR COLORBAR UICONTROL ADDLISTENER
%
% ### Adjust External Colormaps ###
%
% % Example:
% load spine
% image(X)
% brewermap_view({gca})
%
% Very useful! Simply provide a cell array of axes or figure handles when
% calling this function, and their colormaps will be updated in real-time.
% Note that MATLAB versions <=2010 only support axes handles for this.
%
% ### Input and Output Arguments ###
%
% Input Arguments (*==default):
%  N  = NumericScalar, an integer to define the colormap length.
%     = *[], use the maximum defined number nodes for each colorscheme.
%     = {axes/figure handles}, their colormaps will be updated by this function.
%  scheme = StringToken, a ColorBrewer scheme name to select the colorscheme.
%
% Output Arguments (block execution until the figure is deleted!):
%  map    = NumericMatrix, the colormap defined when the figure is closed.
%  scheme = StringToken, the name of the colorscheme given in <map>.
%
% [map,scheme] = brewermap_view(N, scheme)

% ### Input Wrangling ###
%
if nargin<1 || isnumeric(N)&&isempty(N)
	N = 128;
elseif iscell(N)&&numel(N)
	tmp = all(1==cellfun('prodofsize',N)&cellfun(@ishghandle,N));
	assert(tmp,'Input <N> may be a cell array of axes or figure handles.')
else
	assert(isnumeric(N)&&isscalar(N),'Input <N> must be a scalar numeric.')
	assert(isreal(N)&&fix(N)==N&&N>0,'Input <N> must be positive real integer: %g+%gi',N,imag(N))
	N = double(N);
end
%
if nargin<2
	L = brewermap('list');
	scheme = L{1+rem(round(now*1e7),numel(L))};
else
	assert(ischar(scheme)&&isrow(scheme),'Second input <scheme> must be a string.')
end
%
[~,~,h] = bmvUpDt(N, scheme);
%
if nargout
	waitfor(h);
	[~,scheme,map] = bmvUpDt();
end
%
end
%----------------------------------------------------------------------END:colorbrewer_view
function [N,S,ghnd] = bmvUpDt(N,S)
% Draw a new figure or update an existing figure. Callback for buttons & demo.
%
persistent H prvS prvN extH
%
% LHS and RHS slider bounds/limits, and slider step sizes:
lbd = 1;
rbd = 128;
stp = [1,10];
%
[L,V,T] = brewermap('list');
%
switch nargin
	case 0 % Demo initialize OR return colormap
		N = prvN;
		S = prvS;
		if nargout==3 % Return colormap
			ghnd = brewermap(N,H.fRev(S));
			return
		end
	case 1 % Button/Slider callback
		if get(H.demo,'Value') % Demo
			return
		elseif ischar(N) % scheme name
			S = N;
			N = prvN;
		else % parameter N only
			S = prvS;
			N = round(get(H.vSld(1),'Value'));
		end
	case 2 % Demo update OR function call
		if iscell(N) % External axes/figure
			extH = N;
			N = size(colormap(extH{1}),1);
		end
		if isempty(H) || ~ishghandle(H.fig)
			if nargout<3 % Demo
				return
			end
			% Check brewermap version:
			str = 'The function "brewermap" returned an unusual %s.';
			assert(all(35==[numel(L),numel(V),numel(T)]),str,'array size')
			tmp = find(any(diff(+char(T)),2));
			assert(numel(tmp)==2&&all(tmp==[9;17]),str,'scheme name sequence')
			% Create a new figure:
			H = bmvPlot(lbd,rbd,stp,N,S,L);
		else % Update slider position:
			set(H.vSld, 'Value',max(lbd,min(rbd,N)));
			if nargout<3 % Demo
				S = L{1+mod(find(strcmpi(S,L)),numel(L))};
				set(H.bGrp,'SelectedObject',H.bSep(strcmpi(S,L)));
			end
		end
		%
	otherwise
		error('How did this happen? This should not be possible.')
end
%
prvN = N;
prvS = S;
ghnd = H.fig;
%
% Adjust slider if out of range for the selected scheme:
idx = strcmpi(S,L);
idy = strcmp('Qualitative',T{idx}) && V(idx)<N;
if idy
	N = V(idx);
	set(H.vSld, 'Value',max(lbd,min(rbd,N)));
end
%
% Update parameter value text:
set(H.vTxt(1), 'String',sprintf('N = %.0f',N));
%
% Get ColorBrewer colormap and grayscale equivalent:
[map,num,typ] = brewermap(N,H.fRev(S));
mag = sum(map*[0.298936;0.587043;0.114021],2);
%
% Update colorbar values:
set(H.cbAx, 'YLim', [0,abs(N)+(N==0)]+0.5);
set(H.cbIm(1), 'CData',reshape(map,[],1,3))
set(H.cbIm(2), 'CData',repmat(mag,[1,1,3]))
%
% Update 2D line / 3D patch values:
if get(H.D2D3,'Value')
	set(H.ln2D, 'XData',linspace(0,1,abs(N)));
	set(H.ln2D, {'YData'},num2cell([map,mag],1).');
else
	set(H.pt3D, 'XData',map(:,1), 'YData',map(:,2), 'ZData',map(:,3), 'FaceVertexCData',map)
end
%
% Update warning text:
str = {typ;sprintf('%d Nodes',num)};
set(H.warn,'String',str);
%
% Update external axes/figure:
if ~isempty(extH)
	for k = find(cellfun(@ishghandle,extH))
		colormap(extH{k},map);
	end
end
%
end
%----------------------------------------------------------------------END:bmvUpDt
function H = bmvPlot(lbd,rbd,stp,V,S,L)
% Draw a new figure with RGBplot axes, ColorBar axes, and uicontrol sliders.
%
M = 9; % buttons per column
gap = 0.01; % gaps
bth = 0.04; % demo height
btw = 0.09; % demo width
uih = 0.40; % height of UI control group
cbw = 0.24; % width of both colorbars
axh = 1-uih-2*gap; % axes height
wdt = 1-cbw-2*gap; % axes width
%
H.fig = figure('HandleVisibility','callback', 'Color','white',...
	'IntegerHandle','off', 'NumberTitle','off',...
	'Name','ColorBrewer Interactive Scheme Selector');
%
% Add 2D lineplot:
H.ax2D = axes('Parent',H.fig, 'Position',[gap, uih+gap, wdt, axh],...
	'ColorOrder',[1,0,0; 0,1,0; 0,0,1; 0.6,0.6,0.6], 'HitTest','off',...
	'Visible','off', 'XLim',[0,1], 'YLim',[0,1], 'XTick',[], 'YTick',[]);
H.ln2D = line([0,0,0,0;1,1,1,1],[0,0,0,0;1,1,1,1], 'Parent',H.ax2D, 'Visible','off');
%
% Add 3D scatterplot:
H.ax3D = axes('Parent',H.fig, 'OuterPosition',[0, uih, wdt+2*gap, 1-uih],...
	'Visible','on', 'XLim',[0,1], 'YLim',[0,1], 'ZLim',[0,1], 'HitTest','on');
H.pt3D = patch('Parent',H.ax3D, 'XData',[0;1], 'YData',[0;1], 'ZData',[0;1],...
	'Visible','on', 'LineStyle','none', 'FaceColor','none', 'MarkerEdgeColor','none',...
	'Marker','o', 'MarkerFaceColor','flat', 'MarkerSize',10, 'FaceVertexCData',[1,1,0;1,0,1]);
view(H.ax3D,3);
grid(H.ax3D,'on')
xlabel(H.ax3D,'Red')
ylabel(H.ax3D,'Green')
zlabel(H.ax3D,'Blue')
%
% Add warning text:
H.warn = text('Parent',H.ax2D, 'Units','normalized', 'Position',[1,1],...
	'HorizontalAlignment','right', 'VerticalAlignment','top', 'Color','k');
%
% Add demo button:
H.demo = uicontrol(H.fig, 'Style','togglebutton', 'Units','normalized',...
	'Position',[gap,uih+gap+0*bth,btw,bth], 'String','Demo',...
	'Max',1, 'Min',0, 'Callback',@bmvDemo);
% Add 2D/3D button:
H.D2D3 = uicontrol(H.fig, 'Style','togglebutton', 'Units','normalized',...
	'Position',[gap,uih+gap+1*bth,btw,bth], 'String','2D / 3D',...
	'Max',1, 'Min',0, 'Callback',@(h,~)bmv2D3D(H,h));
% Add reverse button:
H.bRev = uicontrol(H.fig, 'Style','togglebutton', 'Units','normalized',...
	'Position',[gap,uih+gap+2*bth,btw,bth], 'String','Reverse',...
	'Max',1, 'Min',0, 'Callback',@(~,~)bmvUpDt());
H.fRev = @(s)[char(42*ones(1,get(H.bRev,'Value'))),s];
%
% Add colorbars:
C = reshape([1,1,1],1,[],3);
H.cbAx(1) = axes('Parent',H.fig, 'Visible','off', 'Units','normalized',...
	'Position',[1-cbw/1,gap,cbw/2-gap,1-2*gap], 'YLim',[0.5,1.5], 'HitTest','off');
H.cbAx(2) = axes('Parent',H.fig, 'Visible','off', 'Units','normalized',...
	'Position',[1-cbw/2,gap,cbw/2-gap,1-2*gap], 'YLim',[0.5,1.5], 'HitTest','off');
H.cbIm(1) = image('Parent',H.cbAx(1), 'CData',C);
H.cbIm(2) = image('Parent',H.cbAx(2), 'CData',C);
%
% Add parameter slider, listener, and corresponding text:
V = max(lbd,min(rbd,V));
H.vTxt = uicontrol(H.fig,'Style','text', 'Units','normalized',...
	'Position',[gap,uih-bth,btw,bth], 'String','X');
H.vSld = uicontrol(H.fig,'Style','slider', 'Units','normalized',...
	'Position',[gap,gap,btw,uih-bth], 'Min',lbd(1), 'Max',rbd(1),...
	'SliderStep',stp(1,:)/(rbd(1)-lbd(1)), 'Value',V(1));
addlistener(H.vSld, 'Value', 'PostSet',@(a,b)bmvUpDt(0));
%
% Add scheme button group:
H.bGrp = uibuttongroup('Parent',H.fig, 'BorderType','none', 'Units','normalized',...
	'BackgroundColor','white', 'Position',[2*gap+btw,gap,wdt-btw-gap,uih-gap]);
% Determine button locations:
Z = 1:numel(L);
Z = Z+(Z>17);
C = (ceil(Z/M)-1)/4;
R = (M-1-mod(Z-1,M))/M;
% Add scheme buttons to group:
for k = numel(L):-1:1
	H.bSep(k) = uicontrol('Style','Toggle', 'String',L{k}, 'Parent',H.bGrp,...
	'Unit','normalized', 'Position',[C(k),R(k),1/4,1/M]);
end
set(H.bGrp,'SelectedObject',H.bSep(strcmpi(S,L)));
set(H.bGrp,'SelectionChangeFcn',@(~,e)bmvUpDt(get(e.NewValue,'String')));
%
end
%----------------------------------------------------------------------END:bmvPlot
function bmv2D3D(H,tgh)
% Switch between 2D-line and 3D-cube: swap visibility and hittest, then update.
%
if get(tgh,'Value')% 2D
	set(H.ax3D, 'HitTest','off', 'Visible','off')
	set(H.ax2D, 'HitTest','on')
	set(H.pt3D, 'Visible','off')
	set(H.ln2D, 'Visible','on')
else % 3D
	set(H.ax2D, 'HitTest','off')
	set(H.ax3D, 'HitTest','on', 'Visible','on')
	set(H.ln2D, 'Visible','off')
	set(H.pt3D, 'Visible','on')
end
bmvUpDt();
%
end
%----------------------------------------------------------------------END:bmv2D3D
function bmvDemo(tgh,~)
% While the toggle button is depressed run a loop showing ColorBrewer schemes.
%
% Initial scheme:
[N,S] = bmvUpDt();
%
% While the toggle button is down, step schemes:
while ishghandle(tgh)&&get(tgh,'Value')
	%
	% Update figure:
	[~,S] = bmvUpDt(round(N),S);
	%
	% Frame-rate faster/slower:
	pause(1.1);
end
%
end
%----------------------------------------------------------------------END:bmvDemo
% Copyright (c) 2015 Stephen Cobeldick
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
% http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and limitations under the License.
%----------------------------------------------------------------------END:license
