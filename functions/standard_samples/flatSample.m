% flatSample.m
%
% Copyright (c) 2018-19, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the
% GNU/GPL-3.0-or-later.
%
% Creates a flat square sample the desired size and distance from the pinhole
% Plate constructs it at the desired level of diffuse scattering. Used either
% by itself or in conjunction with an analytic sphere.
%
% Calling syntax:
%  sample_surface = flatSample(size, dist, diffuse, sigma)
%
% INPUTS:
%  size    - the size of square to make, in mm, e.g. 2mmx2mm
%  dist    - distance from the pinhole plate to the sample, in mm
%  diffuse - level of diffuse scattering, between 0 and 1
%  sigma   - parameter of scattering
%
% OUTPUTS:
%  sample_surface - TriagSirface object containing the surface
function sample_surface = flatSample(size, dist, material)
    V = [-size/2, -dist, -size/2; ...
         -size/2, -dist,  size/2; ...
          0,      -dist,  size/2; ...
          size/2, -dist, -size/2; ...
          size/2, -dist,  size/2];
    F = [1,2,3; ...
         1,3,4; ...
         3,5,4];
    N = [0,1,0; ...
         0,1,0; ...
         0,1,0];

    M = {'default'; 'default'; 'default'};
    mat_lib = containers.Map({'default'}, {material});

    sample_surface = TriagSurface(V, F, N, M, mat_lib);
end

