% switchPlate.m
%
% Copyright (c) 2018-19, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
%
% Switch which pinhole plate model we are using and whether the rays are being
% Generated in Matlab or C then calles the correct gateway function to mex. 
%     Gives enough outputs for most multi-pixel simulations, should be used to
% build up scanning functions. If more outputs are desired, e.g. in single pixel
% simulaitons then the gateway function should be used directly.
%
% Calling syntax:
%  [cnt, killed, numScattersRay] = switch_plate('name', value, ...) 
% 
% INPUTS:
%  plate_represent - How is the pinhole plate being represented
%  sample          - TriagSurface of the sample
%  maxScatter      - The maximum allowed number of scattering events
%  pinhole_surface - TriagSurface of the pinhole plate
%  thePlate        - information on the simple model of the pinhole plate
%  scan_pos        - [scan_pos_x. scan_pos_z]
%  dist            - Distance from the sample to the pinhole plate
%  sphere          - Information on an analytic sphere
%  ray_model       - Are the rays being generated in 'C' or 'MATLAB'
%  which_beam      - What kind of beam is being moddled, 'Effuse, 'Uniform', or
%                    'Gaussian'
%  beam            - Other information on the beam
%
% OUTPUTS:
%  cnt            - The number of detected rays
%  killed         - The number of artificially stopped rays
%  numScattersRay - Histogram of the number of scattering events rays have
%                   undergone before detection
function [cnt, killed, numScattersRay] = switch_plate(varargin)
    
    for i_=1:2:length(varargin)
        switch varargin{i_}
            case 'plate_represent'
                plate_represent = varargin{i_+1};
            case 'sample'
                sample = varargin{i_+1};
            case 'maxScatter'
                maxScatter = varargin{i_+1};
            case 'pinhole_surface'
                pinhole_surface = varargin{i_+1};
            case 'thePlate'
                thePlate = varargin{i_+1};
            case 'scan_pos'
                scan_pos = varargin{i_+1};
            case 'dist'
                dist_to_sample = varargin{i_+1};
            case 'sphere'
                sphere = varargin{i_+1};
            case 'ray_model'
                ray_model = varargin{i_+1};
            case 'which_beam'
                which_beam = varargin{i_+1};
            case 'beam'
                beam = varargin{i_+1};
            otherwise
                warning([' Input ' num2str(i_) ' not recognised.'])
        end
    end
    
    % Switch between the two ways of generating the rays, we either generate in
    % Matlab and get the maximum outputs plus flexibility of source properties
    % or in C for minimal memory useage and slight speed up.
    if strcmp(ray_model, 'MATLAB')
        % We generate the rays in matlab
        
        % TODO: if we have a lot of rays split this up into running ~100,000
        % rays at once and then collating the results.
        
        % Are we doing direct beam or effuse beam?
        if strcmp(which_beam, 'Effuse')
            [ray_pos, ray_dir] = makeEffuse('n_effusive', beam{1}, ...
                'pinhole_c', beam{2}, 'pinhole_r', beam{3}, 'cosine_n', beam{4});
        elseif strcmp(which_beam, 'Uniform') || strcmp(which_beam, 'Gaussian')
            [ray_pos, ray_dir] = makeDirect('n_rays', beam{1}, ...
                'pinhole_c', beam{2}, 'pinhole_r', beam{3}, ...
                'plot_starting', false, 'theta_max', ...
                beam{4}, 'source_model', beam{5}, 'init_angle', beam{6}, ...
                'sigma', beam{7});
        else 
            error('Wrong type of beam to simulate.')
        end
        rays = {ray_pos, ray_dir};

        % We switch which model of the pinhole plate we are using.
        switch plate_represent
            case 'stl'
                [cnt, killed, ~, ~, ~, numScattersRay, ~] = ...
                    traceRays('rays', rays, 'sample', sample, 'maxScatter', ...
                    maxScatter, 'plate', pinhole_surface, 'scan_pos', scan_pos, ...
                    'dist', dist_to_sample, 'sphere', sphere);
            case 'abstract'
                % TODO
            case 'circle'
                [cnt, killed, ~, ~, ~, numScattersRay, ~] = ...
                    traceSimple('rays', rays, 'sample', sample, 'maxScatter', ...
                        maxScatter, 'plate', thePlate, 'scan_pos', scan_pos, ...
                        'dist', dist_to_sample, 'sphere', sphere);
        end
        
        numScattersRay = binMyWay(numScattersRay, maxScatter);
    elseif strcmp(ray_model, 'C')
        % We let C do all the hard work
        switch plate_represent
            case 'stl'
                [cnt, killed, ~, numScattersRay] = traceRaysGen('sample', ...
                    sample, 'maxScatter', maxScatter, 'plate', pinhole_surface, ...
                    'scan_pos', scan_pos, 'dist', dist_to_sample, 'sphere', ...
                    sphere, 'source', which_beam, 'beam', beam);
            case 'abstract'
                % TODO
            case 'circle'
                [cnt, killed, ~, numScattersRay] = traceSimpleGen('sample', ...
                    sample, 'maxScatter', maxScatter, 'plate', thePlate, ...
                    'scan_pos', scan_pos, 'dist', dist_to_sample, 'sphere', ...
                    sphere, 'source', which_beam, 'beam', beam);
        end
    end
end


