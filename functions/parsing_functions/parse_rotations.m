function rot_angles = parse_rotations(inpt)
    str = strtrim(inpt);
    simple_list = ~contains(str, '-');
    cell_lst = strsplit(str(2:end-1), ',');
    if simple_list
        rot_angles = parse_list_input(inpt);
    else
        rot_angles = str2double(cell_lst{1}):str2double(cell_lst{3}):str2double(cell_lst{5});
    end
end