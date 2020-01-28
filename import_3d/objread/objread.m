function [vertices, fdef, fnorm, fmat, materials] = objread(fname)
    % This function parses wavefront object data
    % It reads the mesh vertices, texture coordinates, normal coordinates
    % and face definitions in a .obj file
    %
    % INPUT: fname - wavefront object file full path
    %
    % OUTPUT: vertices - vertex coordinates
    %         fdef - face definitions (index of vertices)
    %         fnorm - normal to each face, in Cartesian [x y z]
    %         fmat - cell array of material key for each face, a string
    %                uniquely identifying a material
    %         materials - material descriptions, as returned by mtlread()
    %
    % Bernard Abayowa, Tec^Edge
    % 11/8/07
    % Modified Dan Seremet, University of Cambridge
    % 21/01/2020

    % set up field types
    vertices = []; textures = []; normals = []; faces = [];
    fmat = {}; material = 'default';

    fid = fopen(fname);

    % parse .obj file
    while true
        tline = fgetl(fid);
        if ~ischar(tline), break; end  % exit at end of file
        ln = sscanf(tline,'%s',1); % extract line type
        switch ln
        case 'mtllib'   % connect with materials library
            mtl_fname = strip(tline(length('mtllib')+1:end));
        case 'o'    % new body definition -- reset current material
            % body = sscanf(tline(2:end), '%s');
            material = 'default';
        case 'usemtl'   % use material for following faces
            material = strip(tline(length('usemtl')+1:end));
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

            fmat{end+1, 1} = material;
            face.v = cellfun(@str2num, data(1, :));
            face.vt = NaN(1, num_vertices); face.vn = NaN(1, num_vertices);
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

    % initialise and construct exported values
    fnorm = zeros(length(faces), 3); fdef = zeros(length(faces), 3);

    for idx = 1:length(faces)
        % fprintf("Face index %d :\t", idx)
        face = faces(idx);
        fdef(idx, :) = face.v;
        fnorm(idx, :) = getNormal(face, normals, vertices);
    end

    % read materials library
    % First must prepend path to filename
    obj_folder = split(string(fname), '/'); obj_folder = obj_folder(1:end-1);
    mtl_fname = fullfile(obj_folder, mtl_fname);
    materials = mtlread(mtl_fname);
end


function normal = getNormal(face, normals, vertices)
    % Determine what the normal to a face is, either by taking
    % the vertex normals from the obj file, if they are specified
    % and are all the same, or by the vector product of two edges.

    % disp(face.v)
    % if all the vertex normals are the same, will take that as the face normal
    if sum(isnan(face.vn), 'all') == 0  % first check if all are present
        ns = normals(face.vn, :);  % take all vertex normals of the face
        if all(ns == ns(1, :))
            disp("Normal taken from data")
            normal = ns(1, :)/norm(ns(1, :));
            return
        end
    end

    % if we got here, means normals info in the file could not be used
    % so calculate it from the vertices. Assumption is that vertices
    % are always in trig order when seen from outside the sample
    v = vertices(face.v, :);
    ab = v(2, :) - v(1, :); ac = v(3, :) - v(1, :);
    normal = cross(ab, ac);
    normal = normal/norm(normal);
end
