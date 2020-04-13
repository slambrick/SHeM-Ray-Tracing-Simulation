% function that plots the outcome of sampling a 2D probability distribution,
% both separately in its two variables, and as a 2D histogram.
% INPUTS : radius - array of values of radius
%          phi - corresponding azimuthal angles
%          rad_limits - range of radius on plots
%          rad_bins - number of bins per side on the 2D plot
%          name - the name of the plot distribution
function plot_distribution_3d(radius, phi, rad_limit, rad_bins, name)
    figure
    pos2d = [radius.*cos(phi), radius.*sin(phi)];
    rad_step = 2*rad_limit / rad_bins;
    rad_edges = -rad_limit:rad_step:rad_limit;

    histogram2(pos2d(:, 1), pos2d(:, 2), rad_edges, rad_edges, 'FaceColor', 'flat',...
    'Normalization', 'probability')

    view([0,0,1]);
    colorbar;
    xlim([-rad_limit, rad_limit])
    ylim([-rad_limit, rad_limit])
    % xlabel([name ' cos(\phi)'])
    % ylabel([name ' sin(\phi)'])
    xlabel('n_{f, x}')
    ylabel('n_{f, y}')
    zlabel('Counts')
    % title([name ' spatial distribution'])
    set(gca, 'FontSize', 18)
    set(gcf, 'Position', [100, 100, 1000, 800])

end