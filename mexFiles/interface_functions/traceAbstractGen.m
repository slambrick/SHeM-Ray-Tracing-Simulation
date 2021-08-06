function [counted, killed, diedNaturally, numScattersRay] = traceAbstractGen(varargin)
    for i_=1:2:length(varargin)
        switch varargin{i_}
            case 'sample'
                sample_surface = varargin{i_+1};
            case 'max_scatter'
                max_scatter = varargin{i_+1};
            case 'plate'
                plate = varargin{i_+1};
            case 'sphere'
                sphere = varargin{i_+1};
            case 'circle'
                circle = varargin{i_+1};
            case 'source'
                which_beam = varargin{i_+1};
            case 'beam'
                beam = varargin{i_+1};
            otherwise
                warning([' Input ' num2str(i_) ' not recognised.'])
        end
    end

    % MATLAB stores matrices by column then row C does row then column. Must
    % take the traspose of the 2D arrays
    % NOTE: it is import these are the right way round
    V = sample_surface.vertices';
    F = int32(sample_surface.faces');
    N = sample_surface.normals';
    B = sample_surface.lattice';
    C = sample_surface.compositions';

    mat_names = sample_surface.materials.keys;
    mat_functions = cell(1, length(mat_names));
    mat_params = cell(1, length(mat_names));
    for idx = 1:length(mat_names)
        mat_functions{idx} = sample_surface.materials(mat_names{idx}).function;
        mat_params{idx} = sample_surface.materials(mat_names{idx}).params;
    end

    % Get the nessacery source information
    switch which_beam
        case 'Uniform'
            source_model = 0;
            theta_max = beam.theta_max;
            sigma_source = 0;
            init_angle = pi*beam.init_angle/180;
        case 'Gaussian'
            source_model = 1;
            theta_max = 0;
            init_angle = pi*beam.init_angle/180;
            sigma_source = beam.sigma_source;
        case 'Effuse'
            source_model = 2;
            theta_max = 0;
            sigma_source = 0;
            init_angle = 0;
        case 'Infinite'
            source_model = 3;
            theta_max = 0;
            sigma_source = 0;
            init_angle = pi*beam.init_angle/180;
    end

    % Pass the source parameters to C
    source_parameters = [beam.pinhole_r, ...
        beam.pinhole_c(1), beam.pinhole_c(2), beam.pinhole_c(3), ...
        theta_max, init_angle, sigma_source];
    
    % Variables passed as structs not objects
    s = sphere.to_struct();
    c = circle.to_struct();
    p = plate.to_struct();
    
    % The calling of the mex function, ... here be dragons ... don't meddle
    % unles you know what you're doing
    [counted, killed, numScattersRay]  = tracingAbstractGenMex(V, F, N, B, C, s, c, p,...
        mat_names, mat_functions, mat_params, max_scatter, beam.n, ...
        source_model, source_parameters);

    % The number of rays that died naturally, rather than being 'killed'
    % because they scattered too many times.
    diedNaturally = beam.n - sum(counted) - killed;
end

