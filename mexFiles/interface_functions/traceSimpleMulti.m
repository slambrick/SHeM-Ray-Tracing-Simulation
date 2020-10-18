% traceSimpleMultiGen.m
%
% Copyright (c) 2020, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
%
% Gatway function for a simple model of the pinhole plate and generating the
% rays in MATLAB.
%
% Calling Syntax:
% [cntr, killed, diedNaturally, numScattersRay] = traceSimpleGen('name', value, ...)
%
% INPUTS:
%  sample     - TriagSurface of the sample
%  maxScatter - The maximum allowed scattering events
%  plate      - Information on the pinhole plate model in a cell array
%  scan_pos   - [scan_pos_x, scan_pos_z]
%  sphere     - Information on the analytic sphere in a cell array
%
% OUTPUTS:
%  cntr           - The number of detected rays
%  killed         - The number of artificailly stopped rays
%  diedNaturally  - The number of rays that did not get detected naturally
%  numScattersRay - Histogram of the number of scattering events detected rays
%                   have undergone
function [cntr, killed, diedNaturally, numScattersRay] = traceSimpleMulti(varargin)
    
    for i_=1:2:length(varargin)
        switch varargin{i_}
            case 'rays'
                ray_pos = varargin{i_+1}{1};
                ray_dir = varargin{i_+1}{2};
            case 'sample'
                sample_surface = varargin{i_+1};
            case 'max_scatter'
                max_scatter = varargin{i_+1};
            case 'plate'
                plate = varargin{i_+1};
            case 'sphere'
                sphere = varargin{i_+1};
            otherwise
                warning([' Input ' num2str(i_) ' not recognised.'])
        end
    end
    
    % MATLAB stores matrices by column then row C does row then column. Must
    % take the traspose of the 2D arrays
    ray_posT = ray_pos';
    ray_dirT = ray_dir';
    VT = sample_surface.vertices';
    FT = int32(sample_surface.faces');
    NT = sample_surface.normals';
    CT = sample_surface.composition;
    
    mat_names = sample_surface.materials.keys;
    mat_functions = cell(1, length(mat_names));
    mat_params = cell(1, length(mat_names));
    for idx = 1:length(mat_names)
        mat_functions{idx} = sample_surface.materials(mat_names{idx}).function;
        mat_params{idx} = sample_surface.materials(mat_names{idx}).params;
    end
    
    s = sphere.to_struct();
    p = plate.to_struct();
    
    % The calling of the mex function, ...
    [cntr, killed, numScattersRay, detected, which_detector]  = ...
        tracingMultiMex(ray_posT, ray_dirT, VT, FT, NT, CT, s, p, mat_names, ...
            mat_functions, mat_params, max_scatter);
    
    % The number of rays that died naturally, rather than being 'killed'
    % because they scattered too many times.
    diedNaturally = size(ray_pos, 1) - cntr - killed;
    
    % Need to remove the excess zeros from these arrays so we don't include the
    % killed rays
    detected = logical(detected);
    which_detector = which_detector(detected);
    numScattersRayDetect = numScattersRay(detected);
    
    % Put the number of scattering events into a histogram for each detector
    numScattersRay = zeros(plate.n_detectors, maxScatter);
    for i_=1:plate.n_detectors
        ind = which_detector == i_;
        numScattersRay(i_,:) = binMyWay(numScattersRayDetect(ind), maxScatter);
    end
    numScattersRay = numScattersRay';   
end

