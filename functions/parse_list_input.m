% parse_list_input.m
%
% Takes a string that is to represent a list of numbers and convert them to 
function lst = parse_list_input(str)
    str = strtrim(str);
    cell_lst = strsplit(str(2:end-1), ',');
    for i_=1:length(cell_lst)
        cell_lst{i_} = str2double(cell_lst{i_});
    end
    lst = cell2mat(cell_lst);
end

