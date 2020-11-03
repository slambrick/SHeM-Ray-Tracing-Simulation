function custom_cmap = custom_colormap(data, fixed_value)
    % Function to generate a custom colourmap where a fixed value
    % is always a fixed color has the same colour regardless
    % of the spread of the data. Courtesy of:
    % https://uk.mathworks.com/matlabcentral/answers/305073-colormap-fixed-middle-value

    L = length(data);

    topColor =    [255, 50, 0] / 256;    % color for maximum data value
    topIntermColor = [255, 200, 0] / 256;
    indexColor =  [1, 1, 1];             % color for fixed data value
    botIntermColor = [0, 0, 255] / 256;
    botColor = [120, 0, 200] / 256;      % color for minimum data value

    % Calculate where proportionally fixed_value lies between minimum and
    % maximum values
    largest = max(data, [], "all");
    smallest = min(data, [], "all");
    index = L*abs(fixed_value-smallest)/(largest-smallest);

    % Create color map ranging from bottom color to index color
    % Multipling number of points by 100 adds more resolution
    customCMap1 = [linspace(botColor(1),botIntermColor(1),50*index)',...
                linspace(botColor(2),botIntermColor(2),50*index)',...
                linspace(botColor(3),botIntermColor(3),50*index)'];

    customCMap2 = [linspace(botIntermColor(1),indexColor(1),50*index)',...
                linspace(botIntermColor(2),indexColor(2),50*index)',...
                linspace(botIntermColor(3),indexColor(3),50*index)'];

    % Create color map ranging from index color to top color
    customCMap3 = [linspace(indexColor(1),topIntermColor(1),50*(L-index))',...
                linspace(indexColor(2),topIntermColor(2),50*(L-index))',...
                linspace(indexColor(3),topIntermColor(3),50*(L-index))'];
    customCMap4 = [linspace(topIntermColor(1),topColor(1),50*(L-index))',...
                linspace(topIntermColor(2),topColor(2),50*(L-index))',...
                linspace(topIntermColor(3),topColor(3),50*(L-index))'];
    custom_cmap = [customCMap1;customCMap2;customCMap3;customCMap4];  % Combine colormaps

end