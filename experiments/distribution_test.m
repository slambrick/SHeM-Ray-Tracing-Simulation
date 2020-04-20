% Copyright (c) 2020, Dan Seremet.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the
% GNU/GPL-3.0-or-later.

function distribution_test()
    clear

    %% compile
    compile_tests()

    %% Parameters
    num_rays = 1e6;
    direction = [sin(pi/6), -cos(pi/6), 0];
    direction = direction/norm(direction);
    normal = [0, 1, 0];

    diff_params = [
        6,   6,   0.1996,... % maxp and maxq
        1,   0,   0,   1,... % rec lattice vectors
        0.0316,  100];       % peak sigma and envelope sigma

    material.function = 'dw_diffraction';
    material.params = [60, 197, 298, 170, 0, diff_params];
    % material.function = 'diffraction';
    % material.params = [0.5, diff_params];
    material.color = [0.8, 0.8, 1.0];

    material

    [theta, phi] = distribution_test_mex(num_rays, direction, material, normal);

    plot_distribution_3d(sin(theta), phi, 1, 100, '\theta')
    plot_distribution_slice(theta, phi, 0, 0.05, 100, '<110> direction diffraction')
end

% figure
% quiver3(-direction(1), -direction(3), -direction(2), ...
%         direction(1), direction(3), direction(2))

% xlabel('X')
% ylabel('Z')
% zlabel('Y')

% origins = repmat([0 0 0], num_rays, 1);
% vectors = [sin(theta).*cos(phi), sin(theta).*sin(phi), cos(theta)];

% hold on
% quiver3(origins(:, 1), origins(:, 2), origins(:, 3), ...
%         vectors(:, 1), vectors(:, 2), vectors(:, 3))
% hold off
