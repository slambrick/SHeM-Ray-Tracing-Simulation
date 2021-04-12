% import_plate.m
%
% Copyright (c) 2018-19, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
%
% Model funciton for importing a new type of pinhole plate. Provide the name of
% the stl CAD file and then adapat the rotation/movement of the pinholep plate
% to put it into the right place.
%
% Calling syntax:
%  pinhole_surface = import_plate(accuracy);
%
% INPUT:
%  plate_fname - name of the stl file including the full path
%
% OUPUT:
%  pinhole_surface - TriagSurface object containing the triangulation of the
%                     pinhole plate.
function pinhole_surface = import_newPlate(accuracy)
    if nargin == 0
        accuracy = 'low';
    end
    
    data = load('objects/defMaterial.mat');
    defMaterial = data.defMaterial;
    materials = containers.Map({'default'}, {defMaterial});
    
    switch accuracy
        case 'low'
            plate_fname = 'pinholePlates/Pinhole_45_PinholePlate_3.0WD_new_forSimulation.stl';
        case 'medium'
            error('not in repository');
        case 'high'
            error('not in repository');
        otherwise
            error('Enter a correct pinhole plate accuracy');
    end
    
    % Import data from file.
    [F, V, N] = stlread(plate_fname);
    fmat = cell(length(F), 1);
    fmat(:) = {'default'};

    % Put inot a TiagSurface object
    flattice = zeros(size(V, 1), 6);
    pinhole_surface = TriagSurface(V, F, N, flattice, fmat, materials);
    
    % Align the plate to match the simulation
    pinhole_surface.plate_align;
end

