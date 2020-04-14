% Copyright (c) 2020, Dan Seremet.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the
% GNU/GPL-3.0-or-later.

% Generate a sample composed of a background and two rows of parallel strips
% made of materials whose properties can be procedurally generated to
% vary with the strip index. The two series of strips are colored as a gradient.

function sample_surface = strip_series(dist_to_sample, working_dist)
    xs = [0, 1];    % x limits of whole sample
    zs = [0, 1];    % z limits of whole sample
    strip_height = 0.01;    % y height of strip above background

    n_strips = 10;  % number of strips on a row
    horiz_fraction = 0.9;  % fraction of x space that a strip occupies
    vert_fraction = 0.8;   % fraction of z space occupied by strips

    x_space = (xs(2) - xs(1)) / n_strips; % x space available for one strip
    x_width = x_space * horiz_fraction; % x space actually occupied by a strip
    z_height = (zs(2) - zs(1)) * vert_fraction / 2;

    avg_zs = (zs(2) + zs(1)) / 2;   % middle of sample z
    row_1_zs = [avg_zs - z_height, avg_zs];   % z limits of first row
    row_2_zs = [avg_zs, avg_zs + z_height];   % z limits of second row

    V = zeros(n_strips*8+2, 3); % vertex coordinates
    F = zeros(n_strips*4+2, 3); % face definitions
    N = repmat([0 1 0], n_strips*4+2, 1);   % face normals
    M = cell(n_strips*4+2, 1);  % face compositions

    % make the background
    [vb, fb] = rectangle(xs, zs, 0);
    V(1:4, :) = vb;
    F(1:2, :) = fb;

    M{1} = 'background';
    M{2} = 'background';

    %% procedurally generate the strips
    for istrip = 1:n_strips
        strip_x_start = (istrip - 1) * x_space;
        strip_xs = [strip_x_start, strip_x_start + x_width];

        vert_idx = (istrip-1)*8 + 4; % index where vertices of current strips start
        face_idx = (istrip-1)*4 + 2; % index where faces of current strips start

        % first row strip
        [v1, f1] = rectangle(strip_xs, row_1_zs, strip_height);

        V((vert_idx + 1):(vert_idx + 4), :) = v1;
        F((face_idx + 1):(face_idx + 2), :) = f1 + vert_idx;
        M{face_idx+1} = ['a' num2str(istrip)];
        M{face_idx+2} = ['a' num2str(istrip)];

        % same for second row strip
        [v2, f2] = rectangle(strip_xs, row_2_zs, strip_height);

        V((vert_idx + 5):(vert_idx + 8), :) = v2;
        F((face_idx + 3):(face_idx + 4), :) = f2 + vert_idx + 4;
        M{face_idx+3} = ['b' num2str(istrip)];
        M{face_idx+4} = ['b' num2str(istrip)];
    end

    %% now create materials library
    matlib = containers.Map();

    mat_bg.function = 'cosine';
    mat_bg.params = [];
    mat_bg.color = [0.4, 0.4, 0.4];
    matlib('background') = mat_bg;

    color_step = 0.8 / n_strips;

    %% procedurally generate materials. Colors are a simple gradient and don't matter
    for istrip = 1:n_strips
        mat_a.function = 'dw_specular';
        mat_a.params = [100 - 8*istrip, 28, 298, 692, 0, 0.1];
        mat_a.color = [0, istrip*color_step, 1 - istrip*color_step];
        matlib(['a' num2str(istrip)]) = mat_a;

        mat_b.function = 'dw_specular';
        mat_b.params = [100 - 8*istrip, 197, 298, 178, 0, 0.1];
        mat_b.color = [0.8, istrip*color_step, 0];
        matlib(['b' num2str(istrip)]) = mat_b;
    end

    %% add arrow to indicate index growing direction and 'sample top'
    V = [V;
         0.3, strip_height, 0.97;
         0.3, strip_height, 0.93;
         0.7, strip_height, 0.95];
    F = [F; length(V) - 2, length(V) - 1, length(V)];
    N = [N; 0, 1, 0];
    M{length(M) + 1} = 'arrow';

    mat_arr.function = 'broad_specular';
    mat_arr.params = [0.5, 0.1];
    mat_arr.color = [1, 0.7, 0.7];

    matlib('arrow') = mat_arr;

    %% create sample object
    sample_surface = TriagSurface(V, F, N, M, matlib);

    % centre the sample first
    sample_surface.moveBy([-(xs(2) + xs(1)) / 2, 0, -avg_zs]);

    % then offset to account for instrument working distance
    sample_surface.moveBy([dist_to_sample - working_dist, - (dist_to_sample + strip_height), 0]);

end

function [vertices, faces] = rectangle(xs, zs, height)
    vertices = [xs(1), height, zs(1);
                xs(2), height, zs(1);
                xs(1), height, zs(2);
                xs(2), height, zs(2);];

    faces = [1, 2, 3;
             4, 3, 2];
end
