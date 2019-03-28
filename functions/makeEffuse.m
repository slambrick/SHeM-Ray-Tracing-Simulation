% Copyright (c) 2018, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
%
% Generates the effuse beam for use in the SHeM simulation. Assumes the effuse
% beam follows a normallty centred cosiune distribution centred on the middle of
% the pinhole.
%
% Calling syntax:
%  [effuse_pos, effuse_dir] = makeEffuse('name', value, ...);
%
% INPUTS:
%  n_effussive - The number of direct beam rays being used
%  pinhole_c   - The location of the centre of the pinhole, 3 element
%  pinhole_r   - The raidus of the pinhole
%  cosine_n    - The exponent of the cosine distribution, cos^nm defaults to 1
%
% OUTPUTS:
%  effuse_pos - 3xn array of the starting positions of the effuse beam rays
%  effuse_dir - 3xn array of the starting directions of the effuse beam ray
function [effuse_pos, effuse_dir] = makeEffuse(varargin)
    
    % Get the inputs
    for i_=1:2:length(varargin)
        switch varargin{i_}
            case 'n_effusive'
                n_effusive = varargin{i_+1};
            case 'pinhole_c'
                pinhole_c = varargin{i_+1};
            case 'pinhole_r'
                pinhole_r = varargin{i_+1};
            case 'cosine_n'
                cosine_n = varargin{i_+1};
            otherwise
                warning([' Input ' num2str(i_) ' not recognised.'])
        end
    end
    
    % Default inputs
    if ~exist('cosine_n', 'var')
        cosine_n = 1;
    end
    
    % Uniformly distributed positions
    pos_phi = 2*pi*rand(n_effusive,1);
    pos_r = pinhole_r*sqrt(rand(n_effusive,1));
    ray_pos2(:,1) = pos_r .* cos(pos_phi);
    ray_pos2(:,3) = pos_r .* sin(pos_phi);
    ray_pos2(:,2) = zeros(n_effusive,1);
    effuse_pos = bsxfun(@plus, ray_pos2, pinhole_c);
    
    % Cos^n distribution of directions
    phis = 2*pi*rand(n_effusive, 1);
    c_theta = (1 - rand(n_effusive, 1)).^(1/(cosine_n + 1));
    s_theta = sqrt(1 - c_theta.^2);
    effuse_dir = [s_theta.*cos(phis), -c_theta, s_theta.*sin(phis)];
end

