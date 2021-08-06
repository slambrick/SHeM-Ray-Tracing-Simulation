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
    elseif strcmp(strtrim(s), 'crude annular')
        % Using the crude annular pinhole plate
        pinhole_model = 'annular';
    elseif strcmp(strtrim(s), 'normal incidence 3mm')
        % Using the normal incidence plate with 3mm disc aperture
        pinhole_model = 'normal';
    elseif strcmp(strtrim(s), 'abstract')
        % An abstract detection direction with a certain size
        pinhole_model = 'abstract';
    elseif strcmp(s(end-3:end), '.stl')
        % Using some ohter pinhole plate model
        pinhole_model = 'stl';
        error('This feature has not yet been implemented');
    else
        % Using a simple pinhole plate model
        pinhole_model = 'N circle';
    end
end

