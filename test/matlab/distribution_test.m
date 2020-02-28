% Copyright (c) 2020, Dan Seremet.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the
% GNU/GPL-3.0-or-later.

clear

%% compile
if true
    compile_tests()
end

%% Parameters
num_rays = 1e5;
direction = [1, -1, 0];
direction = direction/norm(direction);
normal = [0, 1, 0];

material.function = 'broad_specular';
material.params = 0.2;
material.color = [0.8, 0.8, 1.0];

[theta, phi] = distribution_test_mex(num_rays, direction, material, normal);

plot_distribution(theta, phi, [0 pi/2], 50, '\theta')

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
