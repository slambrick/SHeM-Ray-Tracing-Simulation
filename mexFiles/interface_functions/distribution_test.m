% Copyright (c) 2020, Dan Seremet.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the
% GNU/GPL-3.0-or-later.
%
% Interface for testing 
function distribution_test()
    mexCompile(true);

    % Parameters
    num_rays = 1e6;
    direction = [sin(pi/6), -cos(pi/6), 0];
    direction = direction/norm(direction);
    normal = [0, 1, 0];

    diff_params = [
        6,   6,   0.1996,... % maxp and maxq, lambda/a
        cosd(14), sind(14), -sind(14), cosd(14),... % rec lattice vectors
        0.0316,  2];       % peak sigma and envelope sigma

    material.function = 'diffraction';
    material.params = [0.6, diff_params];
    material.color = [0.8, 0.8, 1.0];

    [theta, phi] = distribution_test_mex(num_rays, direction, material, normal);

    plot_distribution_3d(sin(theta), phi, 1, 100, '\theta')
    plot_distribution_slice(theta, phi, 0, 0.09, 200)
end
