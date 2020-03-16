% rotx.m
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
