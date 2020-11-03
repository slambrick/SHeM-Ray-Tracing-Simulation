function sample_surface = photo_stereo_test(working_dist)
    sample_surface = inputSample('fname', 'samples/photoStereoSample.stl', ...
        'working_dist', working_dist, 'sample_dist', working_dist, ...
        'dontMeddle', true);
    
    sample_surface.moveBy([0, -1, 0]);
    sample_surface.moveBy([0, -working_dist, 0])
    sample_surface.rotateY;
end

