function [sample_surface, sphere, circle, sample_description] = sample_import(sample_inputs, sphere, circle, working_dist, dontMeddle, ...
        square_size, init_angle)
    switch sample_inputs.sample_type
        case 'flat'
            circle.make = 0;
            sample_surface = flatSample(square_size, sample_inputs.dist_to_sample, ...
                sample_inputs.material);
            sphere = Sphere(0, sample_inputs.material);
            sample_description = ['A flat square sample size ' ...
                num2str(square_size) 'mm.'];
        case 'strips'
            circle.make = 0;
            sample_surface = strip_series(sample_inputs.dist_to_sample, working_dist);
            sphere = Sphere(0, sample_inputs.material);
            sample_description = 'A sample made up of two series of strips with properties varying across';
        case 'sphere'
            sample_surface = flatSample(square_size, sample_inputs.dist_to_sample, ...
                sample_inputs.material);
            sphere.make = 1;
            circle.make = 0;
            sample_description = ['A single analytic sphere, radius ' ...
                num2str(sphere.radius) 'mm on a flat square of ' num2str(square_size) 'mm.'];
            shift = (working_dist - sample_inputs.dist_to_sample)*tand(init_angle);
            sample_surface.moveBy([-shift, 0, 0]);
            sphere.centre = sphere.centre + [-shift; 0; 0];
        case 'circle'
            sample_surface = flatSample(1, -1000, sample_inputs.material);
            circle.make = 1;
            sphere = Sphere(0, sample_inputs.material);
            sample_description = 'A flat circular sample';
            shift = (working_dist - sample_inputs.dist_to_sample)*tand(init_angle);
            circle.centre = circle.centre + [-shift, 0, 0];
        case 'custom'
            circle.make = 0;
            sample_surface = inputSample('fname', sample_inputs.sample_fname, 'sample_dist', sample_inputs.dist_to_sample, ...
                                         'working_dist', working_dist, 'scale', sample_inputs.scale, ...
                                         'defMaterial', sample_inputs.material, 'dontMeddle', dontMeddle, ...
                                         'init_angle', init_angle);
            sphere = Sphere(0, sample_inputs.material);
            sample_description = sample_inputs.sample_description;
        case 'photoStereo'
            circle.make = 0;
            sample_surface = photo_stereo_test(working_dist);
            sphere.radius = 0.05;
            c = [0.1, -sample_inputs.dist_to_sample - sphere.radius*2/3, 0.1];
            sphere.make = 1;
            sphere = Sphere(1, sample_inputs.material, c, 0.05);
            sample_description = 'Sample containing different features for testsing photometric stereo';
        case 'poly_crystal'
            circle.make = 0;
            sample_surface = inputSample('fname', 'samples/poly_crystal.obj', 'dontMeddle', true);
            sphere = Sphere(0, sample_inputs.material);
            sample_description = 'A polycrystalline sample';
            sample_surface.moveBy([0, -sample_inputs.dist_to_sample, 0]);
        case 'special'
            circle.make = 0;
            % NOTE: I can't remember what this does...
            sample_surface = inputSample('fname', sample_inputs.sample_fname, ...
                'dontMeddle', true, 'scale', 10e-4);
            sample_surface.rotateZ;
            sample_surface.moveBy([0, -2.121, 0]);
            sphere = Sphere(0, sample_inputs.material);
        case 'test_sample'
            circle.make = 0;
            sample_surface = inputSample('fname', 'samples/alternative_design_test_sample.stl', ...
                'working_dist', working_dist', 'scale', 1, 'defMaterial', sample_inputs.material, ...
                'dontMeddle', true, 'init_angle', init_angle);
            sample_surface.rotateGeneral('y', -160); % Rotate this sample
            sample_surface.moveBy([0, -0.2, 0]); % Moves the sample into the y=0 plane
            sphere = Sphere(true, sample_inputs.material, [-0.22, -0.15; ...
                                                0.015,  0.02; ...
                                                0.13, 0.2], [0.03, 0.04]);
            sample_surface.moveBy([tand(init_angle)*(sample_inputs.dist_to_sample - working_dist), -working_dist, 0]);
            sphere = sphere.move([tand(init_angle)*(sample_inputs.dist_to_sample - working_dist); -working_dist; 0]);
            sample_surface.moveBy([0.16, 0, -0.125]);
            sphere = sphere.move([0.16; 0; -0.125]);
            sample_description = 'Two spheres and some blocks';
            
    end

    if dontMeddle
        disp('Sample has not been automatically placed, you will need to do it manually.')
        keyboard
    end

end