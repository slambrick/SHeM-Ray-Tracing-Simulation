function plot_distribution_slice(theta, phi, direction, slice_width, nbins)

    % compute the direction opposite. NB phi is in (-pi, pi)
    % first bring to (0, 2pi) by adding pi,
    % then add another pi and take mod for the opposite
    % then subtract pi to bring back to (-pi, pi)
    % opposite_dir = mod((direction + 2*pi), 2*pi) - pi;

    % convert direction to radians
    direction = direction * pi / 180;

    % filter out the thetas which are within slice_width
    % of the line defined by (cos(direction), sin(direction)) * length
    % flags = 1 for positive direction, -1 for negative direction
    flags = zeros(size(theta));
    for idx = 1:length(theta)
        d_point_line = sin(theta(idx)) * abs(sin(phi(idx) - direction));
        if d_point_line < slice_width
            if cos(phi(idx) - direction) >= 0
                flags(idx) = 1;
            else
                flags(idx) = -1;
            end
        end
    end

    positive_theta = theta(flags == 1);
    negative_theta = - theta(flags == -1);

    figure
    % histogram([negative_theta; positive_theta], nbins, 'BinLimits', [-pi/2 pi/2],...
    %           'Normalization', 'probability', 'EdgeColor', 'None');
    % hold on;

    bin_leap = pi/nbins;
    edges = -pi/2:bin_leap:pi/2;

    [N,edges] = histcounts([negative_theta; positive_theta], edges);
    edges = edges(2:end) - (edges(2)-edges(1))/2;
    plot(edges, N, 'LineWidth', 3);

    % make pretty
    xticks(-pi/2:pi/6:pi/2);
    xticklabels({'-\pi/2', '-\pi/3', '-\pi/6', '0', '\pi/6', '\pi/3', '\pi/2'});
    xlim([-pi/2, pi/2]);

    xlabel('\theta (rad)')
    ylabel('Counts')
end