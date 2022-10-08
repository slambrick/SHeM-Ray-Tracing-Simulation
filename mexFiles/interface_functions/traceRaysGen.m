% Copyright (c) 2018-20, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
%
% Gatway function for a CAD model of the pinhole plate and generating the
% rays in C.
%
% Calling Syntax:
% [cntr, killed, diedNaturally, numScattersRay] = traceRaysGen('name', value, ...)
%
% INPUTS:
%  sample     - TriagSurface of the sample
%  maxScatter - The maximum allowed scattering events
%  plate      - TriagSurface of the pinhole plate
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
function [cntr, killed, diedNaturally, numScattersRay] = traceRaysGen(varargin)
    
    for i_=1:2:length(varargin)
        switch varargin{i_}
            case 'sample'
                sample_surface = varargin{i_+1};
            case 'max_scatter'
                max_scatter = varargin{i_+1};
            case 'plate'
                pinhole_surface = varargin{i_+1};
            case 'sphere'
                sphere = varargin{i_+1};
            case 'circle'
                circle = varargin{i_+1};
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
    FT = int32(sample_surface.faces');
    NT = sample_surface.normals';
    BT = sample_surface.lattice';
    CT = sample_surface.compositions';
    
    VTS = pinhole_surface.vertices';
    FTS = int32(pinhole_surface.faces');
    NTS = pinhole_surface.normals';
    BTS = pinhole_surface.lattice';
    CTS = pinhole_surface.compositions';
    
    mat_names = sample_surface.materials.keys;
    mat_functions = cell(1, length(mat_names));
    mat_params = cell(1, length(mat_names));
    for idx = 1:length(mat_names)
        mat_functions{idx} = sample_surface.materials(mat_names{idx}).function;
        mat_params{idx} = sample_surface.materials(mat_names{idx}).params;
    end
    
    % Need to know how deep the pinhole plate is, how wide it is and how high it
    % is, this is used in determining if rays are detected, this assumes that
    % the pinhole plate is rectangular and is centered on the origin (in x-z)
    backWall = [max(pinhole_surface.vertices(:,2)), ...
        range(pinhole_surface.vertices(:,1)), ...
        range(pinhole_surface.vertices(:,3))];
    
    % Get the nessacery source information
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
    source_parameters = [beam.pinhole_r, ...
        beam.pinhole_c(1), beam.pinhole_c(2), beam.pinhole_c(3), ...
        theta_max, init_angle, sigma_source];
    
    % Variables passed as struct arrays not objects
    [s, n_sphere] = sphere.to_struct();
    c = circle.to_struct();
    
    % The calling of the mex function, ... here be dragons ... don't meddle
    % unles you know what you're doing
    [cntr, killed, numScattersRay]  = ...
        tracingGenMex(VT, FT, NT, BT, CT, VTS, FTS, NTS, BTS, CTS, n_sphere, s, c, backWall, mat_names, ...
                mat_functions, mat_params, max_scatter, beam.n, source_model, source_parameters);
    
    % The number of rays that died naturally, rather than being 'killed'
    % because they scattered too many times.
    diedNaturally = beam.n - cntr - killed;
end

