function addNumericalValues_Figure(data, ampSNR, cycLabel)
    textStrings = num2str(data(:),'%0.2f');  %# Create strings from the matrix values
    textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
    [x,y] = meshgrid(1:8, 1:8);   %# Create x and y coordinates for the strings
    hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
                    'HorizontalAlignment','center');
    midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
    textColors = repmat(data(:) < midValue,1,3);
    textColors(isnan(data(:)),:) = 1;
    set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
    ylabel({'Simulated amplitudes';'(overall SNR)'}); xlabel('Simulated cycles (Simulated abundance)');
    set(gca, 'YTick', 1:8);
    set(gca, 'YTickLabels', ampSNR);
    set(gca, 'XTickLabels', cycLabel);
end