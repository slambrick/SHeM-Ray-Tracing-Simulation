% import_plate.m
%
% Copyright (c) 2018-19, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
%
% Imports the pinhole plate for the Cambridge SHeM and puts it in the right
% position. Only written for that specific pinhole plate, the 3D model of any
% new plate would have to be consulted before it is adapted for that, though
% similar designs will work with only minimal alterations.
%
% Calling syntax:
%  pinhole_surface = import_plate(accuracy);
%
% INPUT:
%  accuracy - there are difference accuracies for the pinhole plate, by default
%             uses the lowest accuracy for speed. Can be 'low', 'medium' or
%             'high'. Only consider not using 'low' if the pinhole plate itself
%             is being investigated.
%
% OUPUT:
%  plate_surface - TriagSurface object containing the triangulation of the
%                  pinhole plate.
function pinhole_surface = import_plate(accuracy)
    if nargin == 0
        accuracy = 'low';
    end
    
    data = load('objects/defMaterial.mat');
    defMaterial = data.defMaterial;
    materials = containers.Map({'default'}, {defMaterial});
    
    switch accuracy
        case 'low'
            plate_fname = 'pinholePlates/pinholePlate_simple1.stl';
        case 'medium'
            plate_fname = 'pinholePlates/pinholePlate_simple2.stl';
        case 'high'
            plate_fname = 'pinholePlates/pinholePlate_simple3.stl';
        otherwise
            error('Enter a correct pinhole plate accuracy');
    end
    
    % Import data from file.
    [F, V, N] = stlread(plate_fname);
    fmat = cell(length(F), 1);
    fmat(:) = {'default'};
    
    % The simulation runs in mm and the object file is in cm
    V = V*10;
    flattice = zeros(size(V, 1), 6);
    % Put inot a TiagSurface object
    pinhole_surface = TriagSurface(V, F, N, flattice, fmat, materials);
    
    % Manipulate into the right position
    pinhole_surface.plate_align;
end

