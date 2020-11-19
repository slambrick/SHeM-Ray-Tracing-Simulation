classdef SimulationInfo < handle
    properties
        counters;
        time;
        time_estimate;
        init_angle;
        beam_param;
        aperture_axes;
        aperture_c;
        dist_to_sample;
        rays_per_pixel;
    end
    
    methods
        function obj = SimulationInfo(time, time_estimate, init_angle, beam_param, ...
                dist_to_sample, n_rays)
            obj.time = time;
            obj.time_estimate = time_estimate;
            obj.init_angle = init_angle;
            obj.beam_param = beam_param;
            obj.aperture_c = NaN;
            obj.aperture_axes = NaN;
            obj.dist_to_sample = dist_to_sample;
            obj.rays_per_pixel = n_rays;
            obj.counters = NaN;
        end
    end
end

