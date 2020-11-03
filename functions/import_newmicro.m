function pinhole_surface = import_newmicro(accuracy)
    if nargin == 0
        accuracy = 'low';
    end
    
    switch accuracy
        case 'low'
            plate_fname = 'pinholePlates/Pinhole_Plate_04_forRayTracing_low.stl';
        case 'medium'
            plate_fname = 'pinholePlates/Pinhole_Plate_04_forRayTracing_medium.stl';
        case 'high'
            plate_fname = 'pinholePlates/Pinhole_Plate_04_forRayTracing_high.stl';
        otherwise
            error('Enter a correct pinhole plate accuracy');
    end
    
    % Import data from file.
    [F, V, N] = stlread(plate_fname);
    
    % Completely diffuse scattering as the 
    C = 1 + zeros(size(F,1), 1);
    P = zeros(size(F,1), 1);
    
    % Put inot a TiagSurface object
    pinhole_surface = TriagSurface(V, N, F, C, P);
    
    % Manipulate into the right position
    pinhole_surface.plate_align_newmicro;
end

