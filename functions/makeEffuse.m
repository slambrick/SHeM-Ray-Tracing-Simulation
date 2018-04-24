function [effuse_pos, effuse_dir, n_effusive] = makeEffuse(n_rays, ...
        effuse_size, pinhole_c)
% Copyright (c) 2018, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
%
% Generates the effuse beam for use in the SHeM simulation. Assumes the effuse
% beam follows a normallty centred cosiune distribution centred on the middle of
% the pinhole.
%
% INPUTS:
%  n_rays      - The number of direct beam rays being used
%  effuse_size - The relative size of the effuse beam to the direct beam, the
%                number of effuse rays is n_rays*effuse_size (rounded)
%  pinhole_c   - The location of the centre of the pinhole
%
% OUTPUTS:
%  effuse_pos - 3xn array of the starting positions of the effuse beam rays
%  effuse_dir - 3xn array of the starting directions of the effuse beam rays
%  n_effusive - The number of effuse rays generated

    n_effusive = round(n_rays*effuse_size);
    
    % The effusive beam is in the centre of the pinhole with a cosine distribution
    effuse_dir = cosineMex(n_effusive, [0 -1 0])';
    effuse_pos = repmat(pinhole_c, n_effusive, 1);
end

