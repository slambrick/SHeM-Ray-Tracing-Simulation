% inputSample.m
%
% Copyright (c) 2018-19, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
%
% Inputs a traingulated surface from a binary .stl file and puts the data into
% a TriagSurface object. The object may then be used as the sample surface in
% the SHeM simulation. 
% 
% By default positions the 'interesting' part of the sample in the middle. The
% interesting part is defined as those vertices that are not the vertices at the
% back and extremes of the sample, i.e. the corners of the block the sample sits
% on. Note that issues may arise if there a holes all the way through the
% sample, in that case the positioning should be done manually.
% 
% Calling syntax:
%  sample_surface = inputSample('name', value, ...)
%
% INPUTS:
%  fname      - string, the name of the .stl file to inport
%  scattering - double, defines the scattering from the sample:
%                1:   diffuse
%                0:   specular
%                2:   broad specular  (requires parameter)
%                0-1: combination of specular and diffuse scattering
%  plate_dist - minimum distance between the sample and the pinhole plate, in
%               mm, defaults to 2.121
%  scale      - if the .stl file is at a different scale to that needed for the
%               simultion. e.g. if the, defaults to 1
%  dontMeddle - bool, if true then don't do any movment/manipulation of the
%               sample, for if this wants to be done manually. Default is 
%               false.
%  parameters - parameter required for the scattering distribution
%
% OUTPUTS:
%  sample_surface - a TriagSurface object containing all the needed information
%                   about the sample
function sample_surface = inputSample(varargin)
    
    for i_=1:2:length(varargin)
        switch varargin{i_}
            case 'fname'
                fname = varargin{i_+1};
            case 'scattering'
                diffuseLvl = varargin{i_+1};
            case 'plate_dist'
                plate_dist = varargin{i_+1};
            case 'scale'
                scale = varargin{i_+1};
            case 'dontMeddle'
                dontMeddle = varargin{i_+1};
            case 'parameters'
                scattering_parameters = varargin{i_+1};
            case 'working_dist'
                working_dist = varargin{i_+1};
            otherwise
                warning([' Input ' num2str(i_) ' not recognised.'])
        end
    end
    
    % Default inputs
    if ~exist('scale', 'var')
        scale = 1;
    end
    
    if ~exist('plate_dist', 'var')
        plate_dist = 2.121;
    end
    
    if ~exist('dontMeddle', 'var')
        dontMeddle = false;
    end
    
    if ~exist('diffuseLvl', 'var')
        diffuseLvl = 1;
    end
    
    if ~exist('scattering_parameters', 'var')
        scattering_parameters = 0;
    end
    
    % Input checking
    if ((diffuseLvl > 1 || diffuseLvl < 0) && diffuseLvl ~= 2)
        error('inputSample: scattering must be between 0 and 1, or be 2.')
    end
    if ~strcmp('.stl', fname(end-3:end))
        error(['inputSample: wrong type of file name provided, only' ...
            '.stl files accepted.'])
    end
    if (plate_dist <= 0)
        error(['inputSample: cannot put the sample on or behind the pinhole' ...
               'plate'])
    end
    
    % Import from file
    [F, V, N] = stlread(fname);
    V = V/scale;
    
    % Set the scattering off the sample and any scattering parameters
    C = diffuseLvl + zeros(size(F,1), 1);
    P = zeros(size(F,1), 1) + scattering_parameters;
    sample_surface = TriagSurface(V, N, F, C, P);
    
    sample_surface.rotateY;
    
    if ~dontMeddle
        % Move the sample into the right starting position
        move = -plate_dist - max(sample_surface.vertices(:,2));
        sample_surface.moveBy([0 move 0]);
        
        % Position the interesting part of the sample into the correct place, so
        % that is is centered at '0 raster movment'
        minY = min(V(:,2));
        backVertices = find(V(:,2) == minY);
        backXs = V(backVertices,1); %#ok<FNDSB>
        ind = ismember(V(:,1), backXs);
        interestingV = V(~ind,:);
        middleX = (max(interestingV(:,1)) + min(interestingV(:,1))) / 2;
        middleZ = (max(interestingV(:,3)) + min(interestingV(:,3))) / 2;
        moveX = (plate_dist - working_dist) - middleX;
        sample_surface.moveBy([moveX, 0, -middleZ]);
    end
end

