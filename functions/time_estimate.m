% Copyright (c) 2018, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.

function t = time_estimate(n_rays, n_effuse, sample_surface, pixels)
% Estimates the time the simulation will take and prints that out the the
% terminal. 
%
% INPUTS:
%  n_rays          - the number of rays per pixel in this simulation
%  sample_surface  - TriagSurface object of the sample
%  pinhole_surface - TriagSurface object of the pinhole plate
%  pixels          - the number of pixels in the simulation
%
% OUTPUT:
%  t - time estimate in seconds
    
    % Time estimate for the direct rays
    t_1 = n_rays*pixels*((sample_surface.nTriag - 3)*1.4962e-07 + 1.8114e-05);
    
    % Time estimate for the effuse rays
    t_2 = n_effuse*pixels*((sample_surface.nTriag - 3)*1.7794e-07 + 1.1557e-05);
    
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
    
    fprintf(['\nThe estimate may be an overestimate when the sample is' ...
             ' a large distance away\nfrom the pinhole plate, or an ' ...
             'underestimate if it is close. The estimation is\nalso likely ' ...
             'to not be accurate for very small simulations or with ' ...
             'non-cosine\nscattering.\n\n']);
    fprintf('Note: this estimate was calculated for an modern laptop quad core i7.\n\n');
    cont = input('Is this okay (y/n)?\n', 's');
    if cont == 'n'
        error('User stopped execution.');
    end
    fprintf('\n');
end

