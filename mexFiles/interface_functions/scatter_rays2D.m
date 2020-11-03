% Copyright (c) 2018, Sam Lambrick.
% All rights reserved.
% This file is part of the Sub-beam Ray Tracing Simulation, subject to the  
% GNU/GPL-3.0-or-later.
%
% Gateway frunction for the C function that scatters rays in 2D. Scatters the
% provided rays off the provided surface returning the final direction positions
% and number of scattering events of the rays. Arrays need to be provided in the
% correct orientation for conversion to C arrays.
%
% INPUTS:
%  surface - struct, {vertices, normals}. Details the 1D surface that the rays
%            are to be scattered off. Both are arrays of two element column
%            vectors, the number of nurmals must be one less than the number of
%            vertices.
%  rays    - struct, {initial_positions, initial_direction}. Details the rays
%            that are to be scattered off the 1D surface. The positions are
%            specified by an array of n 2 element column vectors and the
%            direction is specified by a single 2 element column vector.
%  scattering - struct, {type_of_scattering, scattering_parameters}. Defines the
%               model for the scattering process (specular, diffuse, ...) along
%               with any parameters needed for that scatteirng model. The type
%               of scattering is provided as a string and the parameters as an
%               array.
%
% OUPUTS:
%  final_dir    - array of 2 element column vectors of the final directions of
%                 the arrays after they now longer intersect with the surface.
%  final_pos    - array do 2 element column vectors of the final positions of
%                 the arrays after they no longer intersect with the surface.
%  num_scatters - the number of scattering events that each ray undergoes.
function [final_dir, final_pos, num_scatters] = scatter_rays2D(surface, rays, ...
        scattering)
    % Perform some checks on the inputs
    if size(surface{1}, 1) ~= 2 || size(surface{2}, 1) ~= 2
        error(['The surface element is specified by two element verticies ' ...
            'column vectors and two elemennt surface noraml column vectors.']);
    end
    if size(surface{1}, 2) ~= size(surface{2}, 2) + 1
        error(['The number of surface normals must be one less than the ' ...
            'number of surface vertices.']);
    end
    if size(rays{1}, 1) ~= 2 || size(rays{2}, 1) ~= 2
        error(['The ray direction and positions must be specified by 2 ' ...
            'element column vectors.']);
    end
    
    % The code for the type of scattering
    switch scattering{1}
        case 'specular'
            scattering{1} = 0;
        case 'diffuse'
            scattering{1} = 1;
        case 'uniform'
            scattering{1} = 2;
        case 'broad specular'
            scattering{1} = 3;
        otherwise
            error('Specify an exsisting type of scattering');
    end
    
    % Call the mex function
    [final_dir, final_pos, num_scatters] = scatterRaysMex2D(surface{1}, ...
        surface{2}, rays{1}, rays{2}, scattering{1}, scattering{2});
end

