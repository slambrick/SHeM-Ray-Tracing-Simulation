% Copyright (c) 2018, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.

function sample_surface = inputSample(fname, diffuseLvl, plate_dist, scale, ...
                                      dontMeddle)
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
% INPUTS:
%  fname      - string, the name of the .stl file to inport
%  diffuseLvl - double, a value between 0 and 1 determining how much diffuse
%               scattering there is off of the sample being inported. May also
%               be 2, corresponding to unifomr scattering. Defaults to 
%               completely diffuse, 1.
%  plate_dist - minimum distance between the sample and the pinhole plate, in
%               mm, defaults to 2.121
%  scale      - if the .stl file is at a different scale to that needed for the
%               simultion. e.g. if the, defaults to 1
%  dontMeddle - bool, if true then don't do any movment/manipulation of the
%               sample, for if this wants to be done manually. Default is 
%               false.
%
% OUTPUTS:
%  sample_surface - a TriagSurface object containing all the needed information
%                   about the sample

    if (nargin == 1)
        diffuseLvl = 1;
        plate_dist = 2.121;
        scale = 10;
        dontMeddle = false;
    elseif (nargin == 2)
        plate_dist = 2.121;
        dontMeddle = false;
    elseif (nargin == 3)
        scale = 1;
        dontMeddle = false;
    elseif (nargin == 4) 
        dontMeddle = false;
    end
        
    if ((diffuseLvl > 1 || diffuseLvl < 0) && diffuseLvl ~= 2)
		error('inputSample: diffuseLvl must be between 0 and 1.')
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
    
    % Set the diffuse level
    C = diffuseLvl + zeros(size(F,1), 1);
    sample_surface = TriagSurface(V, N, F, C);
    
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
        moveX = (plate_dist - 2.121) - middleX;
        sample_surface.moveBy([moveX, 0, -middleZ]);
    end
end

