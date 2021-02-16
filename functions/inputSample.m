% inputSample.m
%
% Copyright (c) 2018-20, Sam Lambrick.
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
%  fname      - string, the name of the .stl file to import
%  workingDist- nominal distance to sample for the pinhole design
%  sampleDist - actual minimum distance to the sample surface
%  scale      - if the .stl file is at a different scale to that needed for the
%               simultion. e.g. if the, defaults to 1
%  defMaterial- the default material to be used for stl files or for obj
%               when no material is specified for a face
%  dontMeddle - bool, if true then don't do any movment/manipulation of the
%               sample, for if this wants to be done manually. Default is
%               false.
%
% OUTPUTS:
%  sample_surface - a TriagSurface object containing all the needed information
%                   about the sample
function sample_surface = inputSample(varargin)

    for i_=1:2:length(varargin)
        switch varargin{i_}
            case 'fname'
                fname = varargin{i_+1};
            case 'working_dist'
                plate_dist = varargin{i_+1};
            case 'scale'
                scale = varargin{i_+1};
            case 'dontMeddle'
                dontMeddle = varargin{i_+1};
            case 'sample_dist'
                sample_dist = varargin{i_+1};
            case 'defMaterial'
                defMaterial = varargin{i_+1};
            otherwise
                warning([' Input ' num2str(i_) ' not recognised.'])
        end
    end

    % Default inputs
    if ~exist('scale', 'var')
        scale = 1;
    end
    if ~exist('plate_dist', 'var')
        plate_dist = 1;
    end
    if ~exist('sample_dist', 'var')
        sample_dist = 1;
    end
    if ~exist('dontMeddle', 'var')
        dontMeddle = false;
    end
    if ~exist('defMaterial', 'var')
        data = load('objects/defMaterial.mat');
        defMaterial = data.defMaterial;
    end

    % Import from file
    % Matlab does not allow these 3 trivial operations to be combined :(
    fname_parts = split(fname, '.');
    extension = fname_parts(end);
    extension = extension{1};
    % Why is matlab
    switch extension
        case 'stl'
            [fdef, vertices, fnorm] = stlread(fname);
            fmat = cell(length(fdef), 1); fmat(:) = {'default'};
        case 'obj'
            [vertices, fdef, fnorm, flattice, fmat, materials] = objread(fname);
        otherwise
            error(['File with extension ' extension ' unrecognised.'...
            ' Only stl and obj are supported.'])
    end
    vertices = vertices/scale;

    % add default material to materials set
    if ~exist('materials', 'var')
        materials = containers.Map({'default'}, {defMaterial});
    else
        materials('default') = defMaterial;
    end

    % finally, construct the surface object
    sample_surface = TriagSurface(vertices, fdef, fnorm, flattice, fmat, materials);
    
    if ~dontMeddle
        % Move the sample into the right starting position
        moveY = - sample_dist - max(sample_surface.vertices(:,2));
        sample_surface.moveBy([0 moveY 0]);

        % Position the interesting part of the sample into the correct place, so
        % that is is centered at '0 raster movment'
        V = sample_surface.vertices;
        minY = min(V(:,2));
        backVertices = find(V(:,2) == minY);
        backXs = V(backVertices,1); %#ok<FNDSB>
        ind = ismember(V(:,1), backXs);
        interestingV = V(~ind,:);
        middleX = (max(interestingV(:,1)) + min(interestingV(:,1))) / 2;
        middleZ = (max(interestingV(:,3)) + min(interestingV(:,3))) / 2;
        moveX = (plate_dist - sample_dist) - middleX;
        sample_surface.moveBy([moveX, 0, -middleZ]);
    end
end

