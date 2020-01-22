function obj = readObj(fname)
    %
    % obj = readObj(fname)
    %
    % This function parses wavefront object data
    % It reads the mesh vertices, texture coordinates, normal coordinates
    % and face definitions(grouped by number of vertices) in a .obj file
    %
    %
    % INPUT: fname - wavefront object file full path
    %
    % OUTPUT: obj.v - mesh vertices
    %       : obj.vt - texture coordinates
    %       : obj.vn - normal coordinates
    %       : obj.f - face definition
    %
    % Bernard Abayowa, Tec^Edge
    % 11/8/07
    % Modified Dan Seremet, University of Cambridge
    % 21/01/2020

    % set up field types
    vertices = []; textures = []; normals = []; faces = [];
    material = '';

    fid = fopen(fname);

    % parse .obj file
    while true
        tline = fgetl(fid);
        if ~ischar(tline), break; end  % exit at end of file
        ln = sscanf(tline,'%s',1); % extract line type
        switch ln
            case 'o'    % new body definition -- reset current material
                body = sscanf(tline(2:end), '%s');
                % disp(['Object ', body]);
                material = '';
            case 'usemtl'   % use material for following faces
                material = strip(tline(length('usemtl')+1:end));
                % disp(['Switching material to ', material]);
            case 'v'   % mesh vertices
                vertices = [vertices; sscanf(tline(2:end),'%f')'];
            case 'vt'  % texture coordinate
                textures = [textures; sscanf(tline(3:end),'%f')'];
            case 'vn'  % normal coordinate
                normals = [normals; sscanf(tline(3:end),'%f')'];
            case 'f'   % face definition
                % each face has 3 or more vertex references
                % each vertex reference can have up to 3 'fields' 
                % separated by slashes: vertex/texture/normal.
                str = textscan(tline(2:end),'%s'); str = str{1};
                tokens = split(str, '/');

                num_vertices = length(tokens(:, 1));
                num_fields = length(tokens(1, :));
                data = tokens';

                face.mat = material;
                face.v = cellfun(@str2num, data(1, :));
                face.vt = NaN(num_vertices); face.vn = NaN(num_vertices);
                if num_fields >= 2
                    if ~isempty(data{2, 1})
                        face.vt = cellfun(@str2num, data(2, :));
                    end
                end
                if num_fields >= 3
                    face.vn = cellfun(@str2num, data(3, :));
                end
                faces = [faces; face];
        end
    end
    fclose(fid);

    % set up matlab object
    obj.v = vertices; obj.vt = textures; obj.vn = normals; obj.fs = faces;
