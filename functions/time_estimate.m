% time_estimate.m
%
% Copyright (c) 2018-19, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
%
% Estimates the time the simulation will take and prints that out the the
% terminal.
%
% TODO: update the time estimates. Include the different ray modelling
%
% Calling syntax:
%  t = time_estimate('name', value, ...);
%
% INPUTS:
%  n_rays          - the number of direct beam rays per pixel in this simulation
%  n_effuse        - the number of effuse beam reays per pixel
%  sample_surface  - TriagSurface object of the sample
%  n_pixels        - the number of pixels in the simulation
%  pinhole_model   - the model used for the pinhole, this strongly affects the
%                    time simulations take. 'stl', 'circle', 'abstract'
%
% OUTPUT:
%  t - time estimate in seconds
function t = time_estimate(varargin)
    
    for i_=1:2:length(varargin)
        switch varargin{i_}
            case 'n_rays'
                n_rays = varargin{i_+1};
            case 'n_effuse'
                n_effuse = varargin{i_+1};
            case 'sample_surface'
                sample_surface = varargin{i_+1};
            case 'n_pixels'
                pixels = varargin{i_+1};
            case 'pinhole_model'
                pinhole_model = varargin{i_+1};
            otherwise
                warning([' Input ' num2str(i_) ' not recognised.'])
        end
    end
    
    % Set the constant values for the direct and the effuse rays depending on
    % the pinhole model
    disp(['The pinhole model is: ' pinhole_model]);
    switch pinhole_model
        case 'stl'
            % TODO: update
            ct_direct = 1.4962e-07;
            ct_effuse = 1.7794e-07;
            ct_plate = 1.8114e-05;
        case 'circle'
            % TODO: update
            ct_direct = 1.4962e-07/5;
            ct_effuse = 1.7794e-07/5;
            ct_plate = 0;
        case 'N circle'
            % TODO: update
            ct_direct = 1.4962e-07/5;
            ct_effuse = 1.7794e-07/5;
            ct_plate = 0;
        case 'new micro'
            % TODO: update
            ct_direct = 1.4962e-07;
            ct_effuse = 1.7794e-07;
            ct_plate = 1.8114e-05;
        case 'abstract'
            % TODO
        otherwise
            error(['No time estimate for ' pinhole_model ' pinhole plates.']);
    end
    
    % Time estimate for the direct rays
    t_1 = n_rays*pixels*((sample_surface.nTriag - 3)*ct_direct + ct_plate);
    
    % Time estimate for the effuse rays
    t_2 = n_effuse*pixels*((sample_surface.nTriag - 3)*ct_effuse + ct_plate);
    
    % Total time
    t = (t_1 + t_2)/2;
    
    % If there is only one pixel then only one core can be used
    if pixels == 1
        t = t*4;
    end
    fprintf('\nRough estimate for the time the simulation will take: ~%.0f s\n', t);
    hr = floor(t/(60^2));
    min = floor((t - hr*60*60)/60);
    if min == 60
        min = 0;
        hr = hr + 1;
    end
    fprintf('Which is: ~%i hr %i mins\n\n', hr, min);
    fprintf(['The current time is: ' datestr(now, 'HH:MM') '\n']);
    
    end_hr = mod(str2num(datestr(now, 'HH')) + hr, 24); %#ok<*ST2NM>
    end_min = str2num(datestr(now, 'MM')) + min;
    if end_min >= 60
        x = floor(end_min/60);
        end_min = end_min - x*60;
        end_hr = end_hr + x;
        end_hr = mod(end_hr, 24);
    end
    
    fprintf('So the estimated finish time is: ~%2i:%02i\n', end_hr, end_min);
    
    disp(['The estimate may be an overestimate when the sample is' ...
          ' a large distance away from the pinhole plate, or an ' ...
          'underestimate if it is close. The estimation is also likely ' ...
          'to not be accurate for very small simulations or with ' ...
          'non-cosine\nscattering.\n']);
    disp('Note: this estimate was calculated for an modern laptop quad core i7.\n\n');
    
    % Ask if the user wants to continue if the simulation is more than a couple
    % of minutes
    if (t > 120)
        cont = input('Is this okay (y/n)?', 's');
        if cont == 'n'
            error('User stopped execution.');
        end
    end
    fprintf('\n');
end

