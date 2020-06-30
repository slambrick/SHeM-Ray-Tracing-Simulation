% generate_raster_pattern.m
%
% Copyright (c) 2020, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
%
% Generates the raster pattern for the 2D scan
%
% Calling syntax: 
%  raster_pattern = generate_raster_pattern('name', value, ...)
%
% INPUTS:
%  raster_movement2D - [double, double]m the [x,z] raster movement
%  xrange            - 
%  zrange            - 
%  rot_angle         - 
%  init_angle        - 
%
% OUTPUTS:
%  raster_pattern - 
function raster_pattern = generate_raster_pattern(varargin)
    for i_=1:2:length(varargin)
        switch varargin{i_}
            case 'raster_movment2D'
                raster_movment2D_x = varargin{i_+1}(1);
                raster_movment2D_z = varargin{i_+1}(2);
            case 'xrange'
                xrange = varargin{i_+1};
            case 'zrange'
                zrange = varargin{i_+1};
            case 'rot_angle'
                rot_angle = varargin{i_+1};
            case 'init_angle'
                init_angle = varargin{i_+1};
            otherwise
                warning(['input '  ' not recongnised']);
        end
    end
    
    if ~exist('rot_angle', 'var')
        rot_angle = 0;
        disp('No rotation angle adujstment')
    end
    
    if ~exist('init_angle', 'var')
        init_angle = 0;
    end
    
    if ~exist('init_angle', 'var')
        init_angle = 0;
        disp('No Incident angle specified ignoring incident direction projection.');
    end
    
    % Reduce the xrange appropriatly if there is an incidence angle
    xrange_tmp = xrange*cosd(init_angle);
    
    % Generate the scanning pattern 
    x_pattern = xrange_tmp(1):raster_movment2D_x:xrange_tmp(2);
    z_pattern = zrange(1):raster_movment2D_z:zrange(2);
    [xx, zz] = meshgrid(x_pattern, z_pattern);
    xx = xx(:);
    zz = zz(:);
    
    % Rotate the scanning pattern appropriatly about the y axis to match the
    % sample rotation
    for j_=1:length(xx)
        tmp = roty(rot_angle)*[xx(j_); 0; zz(j_)];
        xx(j_) = tmp(1);
        zz(j_) = tmp(3);
    end
    
    % Explicitly make sure that the mean of the x direction is 0
    xx = xx - mean(xx);
    
    % Stretch the x direction according to correct for the perspective
    xx = xx./cosd(init_angle);
    
    raster_pattern.movement_x = raster_movment2D_x;
    raster_pattern.movement_z = raster_movment2D_z;
    raster_pattern.nx = length(x_pattern);
    raster_pattern.nz = length(z_pattern);
    raster_pattern.xrange = xrange;
    raster_pattern.zrange = zrange;
    raster_pattern.x_pattern = xx;
    raster_pattern.z_pattern = zz;
end

