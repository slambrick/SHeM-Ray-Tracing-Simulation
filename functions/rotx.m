% rotx.m
%
% Copyright (c) 2020, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
%
% Generates a 3D rotation matrix about the x axis
%
% Calling syntax:
%  R = rotx(theta);
%
% INPUTS:
%  theta - angle of rotation (anticlockwise) in degrees
%
% OUTPUT:
%  R - 3x3 rotation matrix
function R = rotx(theta)
    R = [1,           0,            0;
         0, cosd(theta), -sind(theta);
         0, sind(theta),  cosd(theta)];
end
