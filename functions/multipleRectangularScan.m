function datasets = multipleRectangularScan(sample_surface, range_y, raster_movement_y, ...
        xrange, zrange, direct_beam, raster_movement_x, raster_movement_z, ...
        maxScatter, pinhole_surface, effuse_beam, dist_to_sample, ...
        sphere, thePath, pinhole_model, thePlate, apertureAbstract, ray_model, ...
        n_detector)

    % find y positions
    ys = range_y(1):raster_movement_y:range_y(2);
    ny = length(ys);
    datasets = cell(ny, 1);

    % progress bar
    prog = waitbar(0, 'Simulations performed');

    for idx = 1:ny
        y_displacement = ys(idx)

        % move the sample in y by given amount
        surface_copy = copy(sample_surface);
        surface_copy.moveBy([y_displacement, -y_displacement, 0]);
        sphere.c = sphere.c + [y_displacement, -y_displacement, 0];

        % create a directory for this scan
        dirname = sprintf('z%+.4fmm', y_displacement);
        dirname = fullfile(thePath, dirname);
        mkdir(dirname);

        % run the simulation at that z
        datasets{idx} = rectangularScan(surface_copy, xrange, zrange, ...
        direct_beam, raster_movement_x, raster_movement_z, maxScatter, ...
        pinhole_surface, effuse_beam, dist_to_sample, ...
        sphere, dirname, pinhole_model, thePlate, apertureAbstract, ray_model, ...
        n_detector);

        waitbar(idx/ny, prog);
    end

    close(prog); delete(prog);
end