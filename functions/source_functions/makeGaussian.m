% makeDirect.m
%
% Copyright (c) 2018-19, Sam Lambrick & Lubos Vozdecky.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
%
% Create the starting rays according to a beam source mode with a Gaussian 
% virtual source model.
%
% Calling syntax:
%  [ray_pos, ray_dir] = makeEffuse('name', value, ...);
%
% INPUTS:
%  n_rays        - The number of rays being generated
%  pinhole_c     - The loacation of the centre of the pinhole, 3 element
%  pinhole_r     - The radius of the pinhole
%  plot_starting - Should plots of the starting positions of the rays be plotted,
%                  true/false, defaults to false
%  thePath       - Path to save plots, only need if we are plotting
%  init_angle    - The incidence angle, in deg.
%  sigma         - Standard deciation of the Gaussian model of the virtual source
%
% OUTPUTS:
%  ray_pos - 3 column array of rays positions
%  ray_dir - 3 column array of rays directions
function [ray_pos, ray_dir] = makeGaussian(varargin)
    
    for i_=1:2:length(varargin)
        switch varargin{i_}
            case 'n_rays'
                n_rays = varargin{i_+1};
            case 'pinhole_c'
                pinhole_c = varargin{i_+1};
            case 'pinhole_r'
                pinhole_r = varargin{i_+1};
            case 'plot_starting'
                plot_starting_positions = varargin{i_+1};
            case 'thePath'
                thePath = varargin{i_+1};
            case 'init_angle'
                init_angle = varargin{i_+1};
            case 'sigma'
                sigma_source = varargin{i_+1};
            otherwise
                warning([' Input ' num2str(i_) ' not recognised.'])
        end
    end
    
    % Default inputs
    if ~exist('plot_starting_positions', 'var')
        plot_starting_positions = false;
    end
    
    % Perform input checking
    if ~exist('ThePath', 'var') && plot_starting_positions
        error('Cant plot starting positions if no data path provided');
    end
    
    % Uniformly distributed positions
    pos_phi = 2*pi*rand(n_rays,1);
    pos_r = pinhole_r*sqrt(rand(n_rays,1));
    ray_pos2(:,1) = pos_r .* cos(pos_phi);
    ray_pos2(:,3) = pos_r .* sin(pos_phi);
    ray_pos2(:,2) = zeros(n_rays,1);
    
    % The polar angle
    dir_theta = sigma_source*sqrt(-2*log((1-rand(n_rays,1))/1));
    
    % Generate the directions
    dir_phi = 2*pi*rand(n_rays,1);
    ray_dir2 = [cos(dir_theta), sin(dir_theta).*cos(dir_phi), ...
        sin(dir_theta).*sin(dir_phi)];
    
    % Project the ray positions into a the coordinate basis being used for the
    % simulation
    ray_pos = bsxfun(@plus, ray_pos2, pinhole_c);
    
    init_angle = (90 - init_angle)*pi/180;
    ray_dir(:,1) = cos(-init_angle)*ray_dir2(:,1) - sin(-init_angle)*ray_dir2(:,2);
    ray_dir(:,2) = sin(-init_angle)*ray_dir2(:,1) + cos(-init_angle)*ray_dir2(:,2);
    ray_dir(:,3) = ray_dir2(:,3);
    
    % We can plot the starting positions of the rays if we want
    if plot_starting_positions
        plot_ray_positions(ray_pos, thePath)
    end
end

