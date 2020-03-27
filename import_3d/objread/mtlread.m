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
    lineno = 0;

    while true
        lineno = lineno + 1;
        tline = fgetl(fid);
        if ~ischar(tline), break; end  % exit at end of file

        ln = sscanf(tline,'%s',1); % extract line type

        % if the first token is a number, then read into params
        if isempty(ln)
            continue
        else
        if all(ismember(ln, '0123456789+-.eEdD'))
            % if there is no current material, error out
            if ~exist('current', 'var')
                error(['On line ' num2str(lineno) ': Cannot read params line ' tline ' when no current material exists']);
            end
            % cut out any potential comment at the end
            tline = strsplit(tline, '#');
            tline = tline{1};

            % read numbers and add to params
            numbers = sscanf(tline, '%f');
            current.params = [current.params, numbers'];
        else
            % if the first token defines a field, read in that
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
                numbers = sscanf(tline(length('params')+1:end), '%f');
                % transpose to make a row vector
                current.params = numbers';
            end
        end
        end
    end
    materials(name) = current;      % add the last material to the map
end