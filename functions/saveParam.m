function saveParam(thePath, sample_inputs, direct_beam, effuse_beam, ...
        pinhole_plate_inputs, scan_inputs, simulationData)
    paramsFile = 'simulationParameters.txt';
    
    fid = fopen([thePath '/' paramsFile], 'w');
    
    FORMAT1 = '%s = ';
    FORMAT2 = '%2.8f\n';
    FORMAT3 = '%i\n';
    FORMAT4 = '[%f, %f]\n';
    
    fprintf(fid, '%s\n\n', 'Parameters and results from simulation');
    
    fprintf(fid, '%s %s %s\n\n', 'This is a', scan_inputs.type_scan, 'scan.');
    
    fprintf(fid, '%s\n\n', sample_inputs.sample_description);
    
    fprintf(fid, FORMAT1, 'Date');
    fprintf(fid, '%s\n', datestr(now, 'dd-mm-yyyy'));
    
    fprintf(fid, FORMAT1, 'Minimum distance from the pinholePlate to the sample (mm)');
    fprintf(fid, FORMAT2, sample_inputs.dist_to_sample);
    
    fprintf(fid, FORMAT1, 'Design perpendicular working distance (mm)');
    fprintf(fid, FORMAT2, pinhole_plate_inputs.working_dist);
    
    fprintf(fid, FORMAT1, 'Scattering label for the sample [type, parameter]');
    fprintf(fid, FORMAT4, sample_inputs.scattering(1), sample_inputs.scattering(2));
    
    fprintf(fid, FORMAT1, 'Maximum number of allowed sample scatters');
    fprintf(fid, FORMAT3, scan_inputs.maxScatter);
    
    fprintf(fid, FORMAT1, 'Pinhole radius (mm)');
    fprintf(fid, FORMAT2, direct_beam.pinhole_r);
    
    fprintf(fid, FORMAT1, 'Pinhole center (x, z) (mm)');
    fprintf(fid, FORMAT4, direct_beam.pinhole_c(1), direct_beam.pinhole_c(3));
    
    fprintf(fid, FORMAT1, 'Total number of rays per pixel');
    fprintf(fid, FORMAT3, direct_beam.n);
    
    fprintf(fid, FORMAT1, 'The model used for the virtual source');
    fprintf(fid, '%s\n', direct_beam.source_model);
    
    switch direct_beam.source_model
        case 'Uniform'
            fprintf(fid, FORMAT1, 'Maximum angle of the virtual source (rad)');
            fprintf(fid, FORMAT2, direct_beam.theta_max);
        case 'Gaussian'
            fprintf(fid, FORMAT1, 'Standard deviation of the virtaul source (rad)');
            fprintf(fid, FORMAT2, direct_beam.sigma_source);
    end
    
    fprintf(fid, FORMAT1, 'Number of rays in the effusive beam');
    fprintf(fid, FORMAT3, effuse_beam.n);
    
    fprintf(fid, FORMAT1, 'cosine exponant for the effusive beam');
    fprintf(fid, FORMAT3, effuse_beam.cosine_n);
    
    fprintf(fid, FORMAT1, 'Model used for the pinhole plate');
    fprintf(fid, '%s\n', pinhole_plate_inputs.pinhole_model);
    
    fprintf(fid, FORMAT1, 'The number of detectors moddeled');
    fprintf(fid, FORMAT2, pinhole_plate_inputs.n_detectors);
    switch pinhole_plate_inputs.pinhole_model
        case 'stl'
            fprintf(fid, FORMAT1, 'Level of accuracy of pinhole plate');
            fprintf(fid, '%s\n', pinhole_plate_inputs.plate_accuracy);
        case 'circle'
            fprintf(fid, FORMAT1, 'Centre of the detector aperture in the y=0 plane (mm)');
            fprintf(fid, FORMAT4, pinhole_plate_inputs.aperture_c(1), ...
                pinhole_plate_inputs.aperture_c(2));
            fprintf(fid, FORMAT1, 'Axes of the aperture (x,z)/mm');
            fprintf(fid, FORMAT4, pinhole_plate_inputs.aperture_axis(1), ...
                pinhole_plate_inputs.aperture_axis(2));
            fprintf(fid, FORMAT1, 'Is the pinhole plate represented to scatter off?');
            fprintf(fid, FORMAT2, pinhole_plate_inputs.plate_represent);
            fprintf(fid, FORMAT1, 'Diameter of the circular pinhole plate (mm)');
            fprintf(fid, FORMAT2, circle_plate_r);
        case 'N circle'
            fprintf(fid, FORMAT1, 'Centre of the detector apertures in the y=0 plane (mm)');
            fprintf(fid, FORMAT4, pinhole_plate_inputs.aperture_c(1), ...
                pinhole_plate_inputs.aperture_c(2));
            fprintf(fid, FORMAT1, 'Axes of the aperture (x,z)/mm');
            fprintf(fid, FORMAT4, pinhole_plate_inputs.aperture_axis(1), ...
                pinhole_plate_inputs.aperture_axis(2));
            fprintf(fid, FORMAT1, 'Is the pinhole plate represented to scatter off?');
            fprintf(fid, FORMAT2, pinhole_plate_inputs.plate_represent);
            fprintf(fid, FORMAT1, 'Diameter of the circular pinhole plate (mm)');
            fprintf(fid, FORMAT2, circle_plate_r);
        case 'abstract'
            % TODO
    end
    
    switch scan_inputs.type_scan
        case 'rotations'
            fprintf(fid, FORMAT1, 'Angles of rotation');
            fprintf(fid, '[');
            for i_=1:length(scan_inputs.rotationAngles)
                fprintf(fid, '%f,', scan_inputs.rotationAngles(i_));
            end
            
            fprintf(fid, FORMAT1, 'Number of pixels in x');
            fprintf(fid, FORMAT3, simulationData.nx_pixels);

            fprintf(fid, FORMAT1, 'Number of pixels in z');
            fprintf(fid, FORMAT3, simulationData.nz_pixels);

            fprintf(fid, FORMAT1, 'Seperation of pixels in x(mm)');
            fprintf(fid, FORMAT2, scan_inputs.raster_movment2D_x);

            fprintf(fid, FORMAT1, 'Seperation of pixels in z(mm)');
            fprintf(fid, FORMAT2, scan_inputs.raster_movment2D_z);
            
            fprintf(fid, FORMAT1, 'Range of x values (min, max) (mm)');
            fprintf(fid, FORMAT4, scan_inputs.xrange(1), scan_inputs.xrange(2));

            fprintf(fid, FORMAT1, 'Range of z values (min, max) (mm)');
            fprintf(fid, FORMAT4, scan_inputs.zrange(1), scan_inputs.zrange(2));
        case 'rectangular'
            fprintf(fid, FORMAT1, 'Number of pixels in x');
            fprintf(fid, FORMAT3, simulationData.nx_pixels);

            fprintf(fid, FORMAT1, 'Number of pixels in z');
            fprintf(fid, FORMAT3, simulationData.nz_pixels);
            
            fprintf(fid, FORMAT1, 'Seperation of pixels in x(mm)');
            fprintf(fid, FORMAT2, scan_inputs.raster_movment2D_x);

            fprintf(fid, FORMAT1, 'Seperation of pixels in z(mm)');
            fprintf(fid, FORMAT2, scan_inputs.raster_movment2D_z);

            fprintf(fid, FORMAT1, 'Range of x values (min, max) (mm)');
            fprintf(fid, FORMAT4, scan_inputs.xrange(1), scan_inputs.xrange(2));

            fprintf(fid, FORMAT1, 'Range of z values (min, max) (mm)');
            fprintf(fid, FORMAT4, scan_inputs.zrange(1), scan_inputs.zrange(2));
        case 'line'
            fprintf(fid, FORMAT1, 'Seperation of pixels (mm)');
            fprintf(fid, FORMAT2, raster_movment1D);

            fprintf(fid, '%s = %s\n', 'Direction of line scan', Direction);

            fprintf(fid, FORMAT1, 'Range of positions (min, max) (mm)');
            fprintf(fid, '%f, %f\n', range1D(1), range1D(2));
        case 'single pixel'
            fprintf(fid, FORMAT1, 'Total number of pixels');
            fprintf(fid, FORMAT3, simulationData.N_pixels);

            fprintf(fid, FORMAT1, 'Estimate for the time taken (s)');
            fprintf(fid, FORMAT2, simulationData.time_estimate);
    end
    
    fprintf(fid, FORMAT1, 'Time the core simulation took (s)');
    fprintf(fid, FORMAT2, simulationData.time);
    
    fclose(fid);
end

