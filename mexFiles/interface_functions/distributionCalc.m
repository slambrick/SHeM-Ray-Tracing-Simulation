% Copyright (c) 2019, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
%
% Gatway function ...
%
% Calling Syntax:
%
% INPUTS:
%
% OUTPUTS:
function [killed, numScattersRay, final_pos, final_dir] = distributionCalc(varargin)
    for i_=1:2:length(varargin)
        switch varargin{i_}
            case 'sample_surface'
                sample_surface = varargin{i_+1};
            case 'maxScatter'
                maxScatters = varargin{i_+1};
            case 'nrays'
                nrays = varargin{i_+1};
            case 'start_pos'
                start_pos = varargin{i_+1};
            case 'start_dir'
                start_dir = varargin{i_+1};
            otherwise 
                warning([' Input ' num2str(i_) ' not recognised.'])
        end
    end
    
    if size(start_pos, 1) ~= 3
        error('Must be nx3 arrays of positions')
    end
    if size(start_pos) ~= size(start_dir)
        error('Same number of positions and directions must be provided')
    end
    if size(start_pos, 2) ~= nrays && size(start_pos, 2) ~= 1
        error('number of positions must either equal number of rays or 1')
    end
    
    % MATLAB stores matrices by column then row C does row then column. Must
    % take the traspose of the 2D arrays
    VT = sample_surface.vertices';
    FT = int32(sample_surface.faces');
    NT = sample_surface.normals';
    CT = sample_surface.compositions';
    
    mat_names = sample_surface.materials.keys;
    mat_functions = cell(1, length(mat_names));
    mat_params = cell(1, length(mat_names));
    for idx = 1:length(mat_names)
        mat_functions{idx} = sample_surface.materials(mat_names{idx}).function;
        mat_params{idx} = sample_surface.materials(mat_names{idx}).params;
    end
    
    
    [killed, numScattersRay, final_pos, final_dir] = ...
        distributionCalcMex(VT, FT, NT, CT, mat_names, mat_functions, mat_params, ...
                            maxScatters, nrays, start_pos, start_dir);
    
    % Remove the positions and directions of the killed rays
    %ind = numScattersRay == -1;
    %final_pos = final_pos(:, ind);
    %final_dir = final_dir(:, ind);
end

