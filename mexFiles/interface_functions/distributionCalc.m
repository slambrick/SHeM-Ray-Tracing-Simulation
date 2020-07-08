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
    
    % MATLAB stores matrices by column then row C does row then column. Must
    % take the traspose of the 2D arrays
    VT = sample_surface.vertices';
    FT = sample_surface.faces';
    NT = sample_surface.normals';
    CT = sample_surface.composition;
    PT = sample_surface.parameters;
    
    [killed, numScattersRay, final_pos, final_dir] = ...
        distributionCalcMex(VT, FT, NT, CT, PT, maxScatters, nrays, start_pos, start_dir);
    
    % Remove the positions and directions of the killed rays
    %ind = numScattersRay == -1;
    %final_pos = final_pos(:, ind);
    %final_dir = final_dir(:, ind);
end

