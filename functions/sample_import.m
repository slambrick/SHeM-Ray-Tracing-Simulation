function [sample_surface, sphere, sample_description] = sample_import(sample_inputs, sphere, working_dist, dontMeddle)
    switch sample_inputs.sample_type
        case 'flat'
            sample_surface = flatSample(square_size, sample_inputs.dist_to_sample, ...
                sample_inputs.material);
            sphere.make_sphere = 0;
            sample_inputs.sample_description = ['A flat square sample size ' ...
                num2str(square_size) 'mm.'];
        case 'strips'
            sample_surface = strip_series(sample_inputs.dist_to_sample, working_dist);
            sphere.make_sphere = 0;
            sample_description = 'A sample made up of two series of strips with properties varying across';
        case 'sphere'
            sample_surface = flatSample(square_size, sample_inputs.dist_to_sample, ...
                sample_inputs.material);
            sphere.make_sphere = 1;
            sample_description = ['A single analytic sphere, radius ' ...
                num2str(sphere.r) 'mm on a flat square of ' num2str(square_size) 'mm.'];
        case 'custom'
            sample_surface = inputSample('fname', sample_inputs.sample_fname, 'sample_dist', sample_inputs.dist_to_sample, ...
                                         'working_dist', working_dist, 'scale', sample_inputs.scale, ...
                                         'defMaterial', sample_inputs.material, 'dontMeddle', dontMeddle);
            sphere.make_sphere = 0;
            sample_description = sample_inputs.sample_description;
        case 'photoStereo'
            sample_surface = photo_stereo_test(working_dist);
            sphere.make_sphere = 1;
            sphere.r = 0.05;
            sphere.c = [-0.1, -sample_inputs.dist_to_sample - sphere.r*2/3, -0.1];
        case 'special'
            % NOTE: I can't remember what this does...
            sample_surface = inputSample('fname', sample_inputs.sample_fname, ...
                'dontMeddle', true, 'scale', 10e-4);
            sample_surface.rotateZ;
            sample_surface.moveBy([0, -2.121, 0]);
            sphere.make_sphere = 0;
    end

    if dontMeddle
        disp('Sample has not been automatically placed, you will need to do it manually.')
        keyboard
    end

end