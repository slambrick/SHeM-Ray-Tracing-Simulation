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

    material.function = 'diffraction';
    % material.params = [0.9, 0.2];
    material.params = [0.7,...          % diffuse level
                       4, 4, 0.1996,...  % maxp, maxq, ratio
                       1, 0, 0, 1,...   % b1 and b2
                       0.03, 2];        % peak and envelope sigma
    material.color = [0.8, 0.8, 1.0];

    [theta, phi] = distribution_test_mex(num_rays, direction, material, normal);

    plot_distribution_3d(theta, phi, [0 pi/2], 100, '\theta')
    plot_distribution_slice(theta, phi, 0, 0.05, '<110> direction diffraction')
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
