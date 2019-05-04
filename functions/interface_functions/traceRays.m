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
%  dist       - distance from the sample to the pinhole plate
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
            case 'maxScatter'
                maxScatter = varargin{i_+1};
            case 'plate'
                pinhole_surface = varargin{i_+1};
            case 'dist'
                dist_to_sample = varargin{i_+1};
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
    FT = sample_surface.faces';
    NT = sample_surface.normals';
    CT = sample_surface.composition;
    PT = sample_surface.parameters;
    
    VTS = pinhole_surface.vertices';
    FTS = pinhole_surface.faces';
    NTS = pinhole_surface.normals';
    CTS = pinhole_surface.composition;
    PTS = pinhole_surface.parameters;
    
    % Need to know how deep the pinhole plate is, how wide it is and how high it
    % is, this is used in determining if rays are detected, this assumes that
    % the pinhole plate is rectangular and is centered on the origin (in x-z)
    backWall = [max(pinhole_surface.vertices(:,2)), ...
        range(pinhole_surface.vertices(:,1)), ...
        range(pinhole_surface.vertices(:,3))];
    
    % The calling of the mex function, ... here be dragons ... don't fiddle
    [cntr, killed, final_pos, final_dir, numScattersRay, detected]  = ...
        tracingMex(ray_posT, ray_dirT, VT, FT, NT, CT, PT, maxScatter, VTS, FTS, ...
                   NTS, CTS, PTS, backWall, sphere.make, ...
                   sphere.c, sphere.r, sphere.scattering, ...
                   sphere.scattering_parameters);
    
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
end

