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
function pinhole_surface = import_angular
    plate_fname = 'pinholePlates/Pinhole_45_PinholePlate_2.0WD_angularResolution_forSimulation_simpler.stl';
    
    data = load('objects/defMaterial.mat');
    defMaterial = data.defMaterial;
    materials = containers.Map({'default'}, {defMaterial});
    
    % Import data from file.
    [F, V, N] = stlread(plate_fname);
    fmat = cell(length(F), 1);
    fmat(:) = {'default'};
    
    % Put inot a TiagSurface object
    pinhole_surface = TriagSurface(V, F, N, fmat, materials);
    
    % Align the plate to match the simulation
    pinhole_surface.plate_align;
end

