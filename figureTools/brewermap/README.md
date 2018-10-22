BrewerMap
=========

The complete palette of ColorBrewer colormaps for MATLAB. Simple selection by scheme name and map length.


One function provides the complete selection of the ColorBrewer colorschemes, especially intended for mapping and plots with attractive, distinguishable colors.


Simple to use: only the the colormap length and the colorscheme name are needed to select and define an output colormap. The colorscheme can be preselected by the user, after which only the colormap length is required to define an output colormap.

The function can be used as a drop-in replacement for the inbuilt colormap functions and it is compatible with all MATLAB functions that require a colormap. The function consists of just one M-file that provides all of the ColorBrewer colorschemes (no file-clutter!). Downsampling or interpolation of the nodes occurs automatically (if required) using each scheme's complete palette of nodes. As an option, the colormap can be returned reversed.

Calling brewermap('demo') creates a figure that displays all of the ColorBrewer colorschemes.

This product includes color specifications and designs developed by Cynthia Brewer (http://colorbrewer.org/). See the ColorBrewer website for further information about each colorscheme, colorblind suitability, licensing, and citations. This is the only submission on MATLAB File Exchange providing ColorBrewer colorschemes that does NOT break the ColorBrewer license conditions!

### Examples ###

    % Plot a scheme's RGB values:
    rgbplot(brewermap(9,'Blues')) % standard
    rgbplot(brewermap(9,'*Blues')) % reversed
    
    % View information about a colorscheme:
    [~,num,typ] = brewermap(0,'Paired')
    num = 12
    typ = 'Qualitative'
    
    % Multiline plot using matrices:
    N = 6;
    axes('ColorOrder',brewermap(N,'Pastel2'),'NextPlot','replacechildren')
    X = linspace(0,pi*3,1000);
    Y = bsxfun(@(x,n)n*sin(x+2*n*pi/N), X.', 1:N);
    plot(X,Y, 'linewidth',4)
    
    % Multiline plot in a loop:
    N = 6;
    set(0,'DefaultAxesColorOrder',brewermap(N,'Accent'))
    X = linspace(0,pi*3,1000);
    Y = bsxfun(@(x,n)n*sin(x+2*n*pi/N), X.', 1:N);
    for n = 1:N
    plot(X(:),Y(:,n), 'linewidth',4);
    hold all
    end
    
    % New colors for the "colormap" example:
    load spine
    image(X)
    colormap(brewermap([],'*YlGnBu'))
    
    % New colors for the "surf" example:
    [X,Y,Z] = peaks(30);
    surfc(X,Y,Z)
    colormap(brewermap([],'RdYlGn'))
    axis([-3,3,-3,3,-10,5])
    
    % New colors for the "contourcmap" example:
    brewermap('*PuOr'); % preselect the colorscheme.
    load topo
    load coast
    figure
    worldmap(topo, topolegend)
    contourfm(topo, topolegend);
    contourcmap('brewermap', 'Colorbar','on', 'Location','horizontal',...
    'TitleString','Contour Intervals in Meters');
    plotm(lat, long, 'k')

### Note ###

Compared to other functions available on MATLAB File Exchange, this function:
* Consists of just one convenient M-file (no .mat files).
* Requires only the standard ColorBrewer scheme name to select the colorscheme.
* Supports all ColorBrewer colorschemes.
* Uses the complete palette for each scheme (eg. Diverging = 15 colors).
* Outputs the MATLAB standard N-by-3 numeric RGB array.
* Default length is the standard MATLAB default colormap length (same length as the current colormap).
* Is compatible with all MATLAB functions that use colormaps (eg: "contourcmap").
* Includes the option to reverse the colormap color sequence.
* Does not break ColorBrewer's Apache license conditions.

### Note ###

The following files are part of GitHub/git repository, and are not required for using this submission in MATLAB:
* .gitattributes
* README.md
