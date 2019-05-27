function [im, param, beam_param] = formatOutput(simData, dataPath)
    % Core information on the image produced
    [~, im.single{1}] = simData.imageSingle('detector', 2, 'plot', ...
        false);
    [~, im.single{2}] = simData.imageSingle('detector', 2, 'plot', ...
        false);
    im.raster_movement_x = simData.raster_movment_x;
    im.raster_movement_y = simData.raster_movment_z;
    [~, im.multiple{1}] = simData.imageMultiple('detector', 1, 'plot', ...
        false);
    [~, im.mutliple{2}] = simData.imageMultiple('detector', 2, 'plot', ...
        false);
    
    % Main simulation parameters
    param.detector_position{1} = simData.aperture_c(1:2);
    param.detector_position{2} = simData.aperture_c(3:4);
    param.detector_axes{1} = simData.aperture_axes(1:2);
    param.detector_axes{2} = simData.aperture_axes(3:4);
    param.z_sample_to_detector = simData.dist_to_sample;
    param.rays_per_pixel = simData.rays_per_pixel;
    
    % Parameters of the incidence beam
    beam_param.init_angle = simData.init_angle;
    beam_param.source_size = simData.beam_param.theta_max;
    beam_param.source_model = simData.beam_param.source_model;
    beam_param.pinhole_r = simData.beam_param.pinhole_r;
    beam_param.pinhole_c = simData.beam_param.pinhole_c;
    
    % Save the formatted data to a subdirectory
    if ~exist([thePath '/formatted'], 'dir')
        mkdir([dataPath '/formatted'])
    end
    dfile = [dataPath '/formatted/reconstructionSimulation.mat'];
    disp(['Saving to ' dfile]);
    save(dfile, 'im', 'param', 'beam_param')
end

