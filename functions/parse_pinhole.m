function pinhole_model = parse_pinhole(s)
    if strcmpi(strtrim(s), 'cambridge')
        % Using the default cambridge pinhole plate
        pinhole_model = 'stl';
    elseif strcmp(s(end-3:end), '.stl')
        % Using some ohter pinhole plate model
        pinhole_model = 'stl';
        error('This feature has not yet been implemented');
    else
        % Using a simple pinhole plate model
        pinhole_model = 'N circle';
    end
end

