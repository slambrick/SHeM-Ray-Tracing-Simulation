% traceRays.m
%
% Copyright (c) 2018-19, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
%
% Gateway function for a simple model of the pinhole plate having been provided
% with the ray starting positions and directions.
%
% Calling Syntax:
% [cntr, killed, diedNaturally, final_pos, final_dir, ...
%     numScattersRayDetect, numScattersRay] = traceSimple('name', value, ...)
% 
%
% INPUTS:
%  rays       - {ray_pos, ray_dir}
%  sample     - TriagSurface of the sample
%  maxScatter - The maximum allowed scattering events
%  plate      - TraigSurface of the pinhole plate
%  scan_pos   - [scan_pos_x, scan_pos_z]
%  sphere     - Information on the analytic sphere in a cell array
%
%
% OUTPUTS:
%  cntr           - The number of detected rays
%  killed         - The number of artificailly stopped rays
%  diedNaturally  - The number of rays that did not get detected naturally
%  final_pos      - The final positions of all the detected rays
%  final_dir      - the final directions of all the detected rays
%  numScattersRay - The number of scattering events each ray has undergone
%  numScattersRayDetect - The number of scattering events each detected ray
%                         has undergone
function [cntr, killed, diedNaturally, final_pos, final_dir, ...
          numScattersRayDetect, numScattersRay] = traceRays(varargin)
    
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
                pinhole_surface = varargin{i_+1};
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
    BT = sample_surface.lattice';
    CT = sample_surface.compositions';
    
    VTS = pinhole_surface.vertices';
    FTS = int32(pinhole_surface.faces');
    NTS = pinhole_surface.normals';
    BTS = pinhole_surface.lattice';
    CTS = pinhole_surface.compositions';
    
    % Need to know how deep the pinhole plate is, how wide it is and how high it
    % is, this is used in determining if rays are detected, this assumes that
    % the pinhole plate is rectangular and is centered on the origin (in x-z)
    backWall = [max(pinhole_surface.vertices(:,2)), ...
        range(pinhole_surface.vertices(:,1)), ...
        range(pinhole_surface.vertices(:,3))];
    
    mat_names = sample_surface.materials.keys;
    mat_functions = cell(1, length(mat_names));
    mat_params = cell(1, length(mat_names));
    for idx = 1:length(mat_names)
        mat_functions{idx} = sample_surface.materials(mat_names{idx}).function;
        mat_params{idx} = sample_surface.materials(mat_names{idx}).params;
    end
    
    s = sphere.to_struct();
    
    % The calling of the mex function, ...
    [cntr, killed, final_pos, final_dir, numScattersRay, detected]  = ...
        tracingMex(ray_posT, ray_dirT, VT, FT, NT, BT, CT, VTS, FTS, ...
                   NTS, BTS, CTS, s, backWall, mat_names, mat_functions, mat_params, max_scatter);
    
    % The number of rays that died naturally, rather than being 'killed'
    % because they scattered too many times.
    diedNaturally = size(ray_pos, 1) - cntr - killed;
    
    % Need to transpose the results back into the format we want
    final_pos = final_pos';
    final_dir = final_dir';
    
    % Need to remove the excess zeros from these arrays so we don't include the
    % killed rays
    detected = logical(detected);
    final_pos = final_pos(detected,:);
    final_dir = final_dir(detected,:);
    numScattersRayDetect = numScattersRay(detected);
    
    numScattersRayDetect = binMyWay(numScattersRayDetect, max_scatter);
end

