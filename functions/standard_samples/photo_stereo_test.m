function sample_surface = photo_stereo_test(working_dist)
    sample_surface = inputSample('fname', 'simulations/photoStereoSample.stl', ...
        'working_dist', working_dist, 'scattering', 1, 'plate_dist', working_dist, ...
        'dontMeddle', true);
    
    sample_surface.moveBy([0, -1, 0]);
    sample_surface.moveBy([0, -working_dist, 0])
    sample_surface.rotateY;
end

