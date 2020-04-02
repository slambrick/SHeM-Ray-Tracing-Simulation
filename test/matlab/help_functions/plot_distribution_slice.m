function plot_distribution_slice(theta, phi, direction, slice_width, name)

    % compute the direction opposite. NB phi is in (-pi, pi)
    % first bring to (0, 2pi) by adding pi,
    % then add another pi and take mod for the opposite
    % then subtract pi to bring back to (-pi, pi)
    opposite_dir = mod((direction + 2*pi), 2*pi) - pi;

    % filter out the thetas which are within slice_width
    % of the line defined by (cos(direction), sin(direction)) * length
    % flags = 1 for positive direction, -1 for negative direction
    flags = zeros(size(theta));
    for idx = 1:length(theta)
        d_point_line = theta(idx) * abs(sin(phi(idx) - direction));
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
    histogram([negative_theta; positive_theta], 'BinLimits', [-pi/2 pi/2]);
    % title(name)
    set(gca, 'FontSize', 18)
    xlabel('\theta (rad)')
    ylabel('Counts')
end