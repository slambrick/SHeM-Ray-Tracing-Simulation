% traceSimpleGen.m
%
% Copyright (c) 2019, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
%
% Gatway function for a simple model of the pinhole plate and generating the
% rays in C.
%
% Calling Syntax:
% [cntr, killed, diedNaturally, numScattersRay] = traceSimpleGen('name', value, ...)
%
% INPUTS:
%  sample     - TriagSurface of the sample
%  maxScatter - The maximum allowed scattering events
%  plate      - Information on the pinhole plate model in a cell array
%  scan_pos   - [scan_pos_x, scan_pos_z]
%  dist       - distance from the sample to the pinhole plate
%  sphere     - Information on the analytic sphere in a cell array
%  which_beam - What kind of beam is being moddled, 'Effuse, 'Uniform', or
%               'Gaussian'
%  beam       - Information on the beam model in an array
%
% OUTPUTS:
%  cntr           - The number of detected rays
%  killed         - The number of artificailly stopped rays
%  diedNaturally  - The number of rays that did not get detected naturally
%  numScattersRay - Histogram of the number of scattering events detected rays
%                   have undergone
function [cntr, killed, diedNaturally, numScattersRay] = traceSimpleGen(varargin)
    
    for i_=1:2:length(varargin)
        switch varargin{i_}
            case 'sample'
                sample_surface = varargin{i_+1};
            case 'maxScatter'
                maxScatter = varargin{i_+1};
            case 'plate'
                plate = varargin{i_+1};
            case 'dist'
                dist_to_sample = varargin{i_+1};
            case 'sphere'
                sphere = varargin{i_+1};
            case 'source'
                which_beam = varargin{i_+1};
            case 'beam'
                beam = varargin{i_+1};
            otherwise
                warning([' Input ' num2str(i_) ' not recognised.'])
        end
    end
    
    % MATLAB stores matrices by column then row C does row then column. Must
    % take the traspose of the 2D arrays
    VT = sample_surface.vertices';
    FT = sample_surface.faces';
    NT = sample_surface.normals';
    CT = sample_surface.composition;
    PT = sample_surface.parameters;
    
    % Get the nessacery source information
    n_rays = beam.n;
    pinhole_c = beam.pinhole_c;
    pinhole_r = beam.pinhole_r;
    switch which_beam
        case 'Uniform'
            source_model = 0;
            theta_max = beam.theta_max;
            sigma_source = 0;
            init_angle = pi*beam.init_angle/180;
        case 'Gaussian'
            source_model = 1;
            theta_max = 0;
            init_angle = pi*beam.init_angle/180;
            sigma_source = beam.sigma_source;
        case 'Effuse'
            source_model = 2;
            theta_max = 0;
            sigma_source = 0;
            init_angle = 0;
    end
    
    % Pass the source parameters to C
    source_parameters = [pinhole_r, pinhole_c(1), pinhole_c(2), theta_max, ...
        init_angle, sigma_source];
    
    % The calling of the mex function, ... here be dragons ... don't meddle
    [cntr, killed, numScattersRay]  = ...
        tracingSimpleGenMex(VT, FT, NT, CT, PT, maxScatter, ...
                   sphere.make, sphere.c, sphere.r, ...
                   sphere.scattering, sphere.scattering_parameter, plate{1}, ...
                   plate{3}, plate{4}, ...
                   plate{5}, n_rays, source_model, source_parameters);
    
    
    % The number of rays that died naturally, rather than being 'killed'
    % because they scattered too many times.
    diedNaturally = n_rays - cntr - killed;
end

