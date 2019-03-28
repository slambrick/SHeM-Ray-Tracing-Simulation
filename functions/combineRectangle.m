% combineRectangle.m
%
% Copyright (c) 2018-19, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
%
% Takes two RectangleInfo objects and combines them. The two simulations
% need to be identical in every respect (the function will test this, but the
% tests are not exhaustive), with the exception of the number of pixels. Used
% for lengthy simulations parts of which are performed on different computers or
% on different nodes in a cluster.
%     For lengthy simulations it is recommended to split the simulation into
% multiple simulations each with fewer rays per pixel and then combine them
% to gain the desired signal-to-noise. This will allow the simulation to be run
% on multiple machines at the same time, or multiple nodes on a cluster at the
% same time.
%
% Calling syntax:
%  outData = combineRectangle(data1, data2)
%
% INPUTS:
%  data1 - RectangleInfo of one simulation
%  data2 - RectangleInfo of another simulation
%
% OUTPUT:
%  outData - RectangleInfo object of the combined data
function outData = combineRectangle(data1, data2)

    % The two simulations need to have the same number of pixels
    if (data1.nx_pixels ~= data2.nx_pixels) || (data1.nz_pixels ~= data2.nz_pixels)
        error('Two data sets not compatible.')
    end
    
    % The two simulations need to have the same range of x & z values
    % Safe to compare these floats because they have undergone no floating
    % point operations.
    if ~isequal(data1.xrange, data2.xrange) || ~isequal(data1.zrange, data2.zrange)
        error('Two data sets not compatible')
    end
    
    % The simulation data contains a sample surface, check that they are
    % the same. Be careful because of comparing floats.
    fprintf(['Double check that the two data sets were of exactly the' ...
        ' same sample in the same conditions.\n'])
    fprintf('The data sets contain the triangulated object. Plot and check.\n')
    keyboard
    
    n_rays = data1.rays_per_pixel + data2.rays_per_pixel;
    time = data1.time + data2.time;
    time_estimate = data1.time_estimate + data2.time_estimate;
    counters = data1.counters + data2.counters;
    num_killed = data1.num_killed + data2.num_killed;
    cntr_effuse = data1.counter_effusive + data2.counter_effusive;
    
    outData = RectangleInfo(counters, num_killed, data1.sample_surface, ...
        data1.xrange, data1.zrange, data1.raster_movment, n_rays, time, ...
        time_estimate, cntr_effuse);
end

