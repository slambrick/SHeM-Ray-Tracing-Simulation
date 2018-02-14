% Copyright (c) 2018, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.

function ray_dir = random_dir(n_rays, alpha)
% random_dir2.m
%
% Generates random directions for the rays using p.d.f. of alpha*x^2. For use in
% the simulation the coordinates must be transformed so that the rays are
% trvelling out of the pinhole, this function centres the directions aroun
% [0,-1,0].
% 
% INPUTS
%  n_rays - the number of directions to generate
%  alpha  - the constant in the p.d.f. (optional)

    phi = 2*pi*rand(n_rays, 1);
    
    alpha = alpha*1000^3;
    
    if nargin == 1
        % in rad^-2
        alpha = 69715282455;
    end
    
    R = rand(n_rays, 1);
    
    theta = (3*R./repmat(alpha,n_rays,1)).^(1/3);
    
    % Generate the direction according to the desired distribution
    extra_dir = [sin(theta).*cos(phi), zeros(n_rays, 1), sin(theta).*sin(phi)];
    
    % Need to add extra_dir to the main direction
    dir = repmat([0 -1 0], n_rays, 1);
    ray_dir = dir + extra_dir;
end

