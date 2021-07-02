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
%  max_scatter      - The maximum allowed number of scattering events
%  pinhole_surface - TriagSurface of the pinhole plate
%  thePlate        - information on the simple model of the pinhole plate
%  scan_pos        - [scan_pos_x. scan_pos_z]
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
            case 'max_scatter'
                max_scatter = varargin{i_+1};
            case 'pinhole_surface'
                pinhole_surface = varargin{i_+1};
            case 'thePlate'
                thePlate = varargin{i_+1};
            case 'sphere'
                sphere = varargin{i_+1};
            case 'circle'
                circle = varargin{i_+1};
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
        
        % I see no need to keep this part of the code, can abandon it and
        % keep C ray generation only
        
        % TODO: if we have a lot of rays split this up into running ~100,000
        % rays at once and then collating the results.
        
        % Which model of the beam are we using
        switch which_beam
            case 'Effuse'
                [ray_pos, ray_dir] = makeEffuse('n_effusive', beam.n, ...
                    'pinhole_c', beam.pinhole_c, 'pinhole_r', beam.pinhole_r, ...
                    'cosine_n', beam.cosine_n);
            case 'Uniform'
                [ray_pos, ray_dir] = makeUniform('n_rays', beam.n, ...
                    'pinhole_c', beam.pinhole_c, 'pinhole_r', beam.pinhole_r, ...
                    'plot_starting', false, 'theta_max', ...
                    beam.theta_max, 'init_angle', beam.init_angle);
            case 'Gaussian'
                [ray_pos, ray_dir] = makeGaussian('n_rays', beam.n, ...
                    'pinhole_c', beam.pinhole_c, 'pinhole_r', beam.pinhole_r, ...
                    'plot_starting', false, 'init_angle', beam.init_angle, ...
                    'sigma', beam.sigma_source);
            otherwise 
                error('Wrong type of beam to simulate.')
        end
        rays = {ray_pos, ray_dir};

        % We switch which model of the pinhole plate we are using.
        switch plate_represent
            case 'stl'
                [cnt, killed, ~, ~, ~, numScattersRay, ~] = ...
                    traceRays('rays', rays, 'sample', sample, 'max_scatter', ...
                    max_scatter, 'plate', pinhole_surface, ...
                    'sphere', sphere, 'circle', circle');
            case 'N circle'
                [cnt, killed, ~, numScattersRay] = traceSimpleMulti('rays', rays, ...
                    'sample', sample', 'max_scatter', max_scatter, ...
                    'sphere', sphere, 'circle', circle', 'plate', thePlate);
            case 'abstract'
                % TODO
        end
        
    elseif strcmp(ray_model, 'C')
        % We let C do all the hard work
        switch plate_represent
            case 'stl'
                [cnt, killed, ~, numScattersRay] = traceRaysGen('sample', ...
                    sample, 'max_scatter', max_scatter, 'plate', pinhole_surface, ...
                    'sphere', sphere, 'circle', circle', 'source', which_beam, 'beam', beam);
            case 'abstract'
                % TODO
            case 'N circle'
                [cnt, killed, ~, numScattersRay] = traceSimpleMultiGen('sample', ...
                    sample, 'max_scatter', max_scatter, 'plate', thePlate, ...
                    'sphere', sphere, 'circle', circle', 'source', which_beam, 'beam', beam);
        end
    end
end


