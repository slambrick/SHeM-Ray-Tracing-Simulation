% Copyright (c) 2018, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.

function sample_surface = composition_sample(dist_to_sample)
% This function constructs a flat sample with four regions with different
% ratios of specular to diffuse scattering. Returns a TriagSurface object. The
% four regions of the sample have 1%, 2%, 3%, 4% specular scattering while the
% background has completely diffuse scattering. The regions are 400um square,
% seperated by 200um. 
%     It should be straight forward to modify this function to test the
% detection from different scattering distributions.
%
% The single argument is the sample-pinhole plate distance.
    if nargin == 0
        dist_to_sample = 2.121;
    end
    
    V = [  4, 0,    4;
           4, 0,   -4;
         0.5, 0,  0.5;
         0.5, 0,  0.1;
         0.5, 0, -0.1;
         0.5, 0, -0.5;
         0.1, 0,  0.5;
         0.1, 0,  0.1;
         0.1, 0, -0.1;
         0.1, 0, -0.5;
        -0.1, 0,  0.5;
        -0.1, 0,  0.1;
        -0.1, 0, -0.1;
        -0.1, 0, -0.5;
        -0.5, 0,  0.5;
        -0.5, 0,  0.1;
        -0.5, 0, -0.1;
        -0.5, 0, -0.5;
          -4, 0,    4;
          -4, 0,   -4];
    
    F = [1,  2,  4;
         1,  3,  4;
         1,  7,  3;
         2,  4,  5;
         2,  5,  6;
         2,  6, 10;
         2, 10, 20;
         1,  7, 19;
         7, 11, 19;
        11, 15, 19;
        15, 16, 19;
        16, 19, 20;
        16, 17, 20;
        17, 18, 20;
        14, 18, 20;
        10, 14, 20;
         4,  5,  9;
         4,  8,  9;
         8,  9, 13;
         8, 12, 13;
         7,  8, 12;
         7, 11, 12;
         9, 10, 14;
         9, 13, 14;
        12, 13, 17;
        12, 16, 17;
         3,  4,  8;
         3,  7,  8;
         5,  6, 10;
         5,  9, 10;
        11, 12, 16;
        11, 15, 16;
        13, 14, 18
        13, 17, 18];
    
    N = repmat([0 1 0], 34, 1);
    
    C1 = repmat([0.98], 26, 1); %#ok<NBRAK>
    C = [C1; 0.99; 0.99; 1; 1; 0.97; 0.97; 0.96; 0.96];
    
    sample_surface = TriagSurface(V, N, F, C);
    
    sample_surface.moveBy([dist_to_sample - 2.121, -dist_to_sample, 0]);
end

