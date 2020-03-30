% formatOutputRotation.m
%
% Copyright (c) 2019-20, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
%
% Formats the results of a simulation into a more useful/simple format when
% there are a series of rotated simulations.
% 
% Calling syntax:
%  [im, param, beam_param] = formatOutputRotations(simData, dataPath, rot_angles)
% 
% INPUTS:
%  simData    - Cell array of simulation data objects
%  dataPath   - Relative path to save the formatted output to. Output is
%               put in a subdirectory in this path
%  rot_angles - Vector of the rotation angles used.
%
% OUTPUTS:
%  im         - Cell array of image results
%  param      - Cell array of general parameters of the simulation
%  beam_param - Cell array of parameters for the set up of the beam
function [im, param, beam_param] = formatOutputRotation(simData, dataPath, rot_angles)
    for i_=1:length(simData)
        save_path = [dataPath '/rotation' num2str(rot_angles(i_))];
        if ~exist(save_path, 'dir')
            mkdir(save_path)
        end
        [im{i_}, param{i_}, beam_param{i_}] = ...
            simData{i_}.formatOutput(save_path); %#ok<AGROW>
    end
end

