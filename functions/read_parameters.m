% read_parameters.m
%
% 
%
% Reads the parameter file and provides a cell array of the read parameters.
% Does not convert the types or assign to variables.
function param_list = read_parameters(param_fname)

    fid = fopen(param_fname);

    param_list = {};

    tline = fgetl(fid);
    i_ = 1;
    while ischar(tline)
        if ~isempty(tline)
            if tline(1) == '%'
                % Line is a comment, ignore
            else
                % Line is not a comment, get whatever comes after the colon
                split_line = strsplit(tline, ':');
                param_list{i_} = split_line{2}; %#ok<AGROW>
                i_ = i_ + 1;
            end
        end
        tline = fgetl(fid);
    end

    fclose(fid);

end

