function materials = mtlread(mtl_fname)
% Read the data stored in a .mtl file associated with an .obj file.
% This data describes the materials referred to from the obj file.
% This is not a general mtl parser, it only extracts properties we care about:
% name, (diffuse) colour, scattering function and parameters.
%
% INPUT: mtl_fname - the filename of the materials library
%
% OUTPUT: a containers.Map where the keys are:
%   - name (string), a unique key
%         and the values are structs with:
%   - function (string), a key to a scattering function
%   - params (double vector), the parameters appropriate for the given function
%   - color (3x1 double) the material colour in [r g b]
%
    materials = containers.Map('KeyType', 'char', 'ValueType', 'any');
    fid = fopen(mtl_fname);

    while true
        tline = fgetl(fid);
        if ~ischar(tline), break; end  % exit at end of file
        ln = sscanf(tline,'%s',1); % extract line type
        switch ln
        case 'newmtl'   % add a new material
            % put the previously read material into the map
            if exist('name', 'var') && exist('current', 'var')
                materials(name) = current;
            end
            % and reset name and current material
            name = strip(tline(length('newmtl')+1:end));
            current = struct('function', '', 'params', [], 'color', [1, 1, 1]);
        case 'Kd'
            current.color = sscanf(tline(3:end), '%f');
            current.color = current.color';
        case 'func'
            current.function = strip(tline(length('func')+1:end));
        case 'params'
            current.params = sscanf(tline(length('params')+1:end), '%f');
        end
    end
    materials(name) = current;      % add the last material to the map
end