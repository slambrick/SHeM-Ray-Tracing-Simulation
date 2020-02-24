% function that plots the outcome of sampling a 2D probability distribution,
% both separately in its two variables, and as a 2D histogram.
% INPUTS : radius - array of values of radius
%          phi - corresponding azimuthal angles
%          rad_limits - range of radius on plots
%          rad_bins - number of bins per side on the 2D plot
%          name - the name of the plot distribution
function plot_distribution(radius, phi, rad_limits, rad_bins, name)
    figure
    histogram(phi, 'BinLimits', [-pi, pi])
    title('\phi distribution')

    figure
    histogram(radius, 'BinLimits', rad_limits)
    title([name ' distribution'])

    figure
    pos2d = [radius.*cos(phi), radius.*sin(phi)];
    max_rad = max(rad_limits);
    rad_step = 2*max_rad / rad_bins;
    rad_edges = -max_rad:rad_step:max_rad;

    hist3(pos2d, 'Edges', {rad_edges rad_edges}, 'FaceColor','interp','CDataMode','auto')
    xlabel('X')
    ylabel('Z')
    title([name ' spatial distribution'])

end