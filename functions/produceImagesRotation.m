% produceImagesRotation.m
%
% Copyright (c) 2020, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
%
% Loop through a whole list of rotated simulations to produce the images. Does
% not compensate for the scan pattern changes.
%
% Calling syntax;
%  produceImagesRotation(simulationData, thePath, rot_angles)
%
% INPUTS:
%  simulationData - cell array of rectangleInfo objects
%  thePath        - data directory to save the images to
%  rot_angles     - list of the rotation angles used
function produceImagesRotation(simulationData, thePath, rot_angles)
    if length(simulationData) ~= length(rot_angles)
        error('Number of rotations must match the number of data sets');
    end
    
    for i_=1:length(simulationData)
        subPath = [thePath '/rotation' num2str(rot_angles(i_))];
        if ~exist(subPath, 'dir')
            mkdir(subPath)
        end
        
        simulationData{i_}.produceImages(subPath);
    end
end

