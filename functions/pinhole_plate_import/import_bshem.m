% import_plate.m
%
% Copyright (c) 2025, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
%
% Imports the pinhole plate for the Cambridge B-SHeM and puts it in the right
% position. Only written for that specific pinhole plate, the 3D model of any
% new plate would have to be consulted before it is adapted for that, though
% similar designs will work with only minimal alterations.
%
% Calling syntax:
%  pinhole_surface = import_plate();
%
% OUPUT:
%  plate_surface - TriagSurface object containing the triangulation of the
%                  pinhole plate.
function pinhole_surface = import_bshem()
    data = load('objects/defMaterial.mat');
    defMaterial = data.defMaterial;
    materials = containers.Map({'default'}, {defMaterial});
   
    plate_fname = 'pinholePlates/bshem_pinholeplate_diffraction_simulation2.stl';
    
    % Import data from file.
    [F, V, N] = stlread(plate_fname);
    fmat = cell(length(F), 1);
    fmat(:) = {'default'};
    
    % The simulation runs in mm and the object file is in cm
    V = V*10;
    flattice = zeros(size(V, 1), 6);
    % Put inot a TiagSurface object
    pinhole_surface = TriagSurface(V, F, N, flattice, fmat, materials);
    %  pinhole_surface.rotateX(270);
    pinhole_surface.moveBy([-2.25, -2+0.00667807, 0]);

    % Manipulate into the right position
    %pinhole_surface.plate_align;
end

