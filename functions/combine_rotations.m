% combine_rotations.m
%
% Copyright (c) 2020, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
%
% Combines two cell array lists of rotations simulations preformed.
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
function outData = combine_rotations(simData1, simData2)
    % TODO: add in checks on parameter structs to ensure the simulations are
    % compatibel
    
    if length(simData1) ~= length(simData2)
        error('Different number of rotations in each data file.')
    end

    outData = cell(size(simData1));
    for i_=1:length(simData1)
        outData{i_} = combineRectangle(simData1{i_}, simData2{i_});
    end
end

