% rotz.m
%
% Generates a 3D rotation matrix about the z axis
%
% Calling syntax:
%  R = rotz(theta);
%
% INPUTS:
%  theta - angle of rotation (anticlockwise) in degrees
%
% OUTPUT:
%  R - 3x3 rotation matrix
function R = rotz(theta)
    R = [cosd(theta), -sind(theta), 0;
         sind(theta),  cosd(theta), 0;
                   0,            0, 1]
end
