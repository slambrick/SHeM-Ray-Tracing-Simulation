function pinhole_model = parse_pinhole(s)
    if strcmpi(strtrim(s), 'cambridge')
        % Using the default cambridge pinhole plate (old chamber)
        pinhole_model = 'stl';
    elseif strcmpi(strtrim(s), 'cambridge new')
        % New chamber pinhole plate
        pinhole_model = 'new';
    elseif strcmpi(strtrim(s), 'angular resolution')
        % Using the angular resolution pinhole plate
        pinhole_model = 'angular';
    elseif strcmp(s(end-3:end), '.stl')
        % Using some ohter pinhole plate model
        pinhole_model = 'stl';
        error('This feature has not yet been implemented');
    else
        % Using a simple pinhole plate model
        pinhole_model = 'N circle';
    end
end

