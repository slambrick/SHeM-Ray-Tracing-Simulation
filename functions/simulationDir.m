% simulationDir.m
%
% Copyright (c) 2018-19, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
%
% Calling syntax:
%  thePath = simulationDir(name)
%
% INPUT:
%  name - string label of the simulation folder
%
% OUTPUT:
%  thePath - full simulation path to store simulation results in 
function thePath = simulationDir(name)
    % Get all th esubdirectories in the simulations folder
    d = dir('simulations');
    isub = [d(:).isdir]; %# returns logical vector
    nameFolds = {d(isub).name}';
    nameFolds(ismember(nameFolds,{'.','..'})) = [];
    
    % Get the largest index of folders so far in the simulations directory
    dir_n = max_n(nameFolds);
    dir_n = dir_n + 1;
    
    % Create the name of the directory to save simulation results to
    dir_n = sprintf('%04i', dir_n);
    thePath = [dir_n '_' name];
end

% Get the maximum index of the cell array of folders.
function n = max_n(nameFolds)
    if length(nameFolds) == 0
        n = 1;
    else
        n = 0;
        for i_=1:length(nameFolds)
            this_n = str2num(nameFolds{i_}(1:4));
            if this_n > n
                n = this_n;
            end
        end
    end
end
