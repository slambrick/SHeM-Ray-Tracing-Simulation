% parse_yes_no.m
%
% Takes a string and returns true if the string is 'yes' or 'true' and false if
% it is 'no' or 'false' in any mixture of upper and lower case letters and with
% any leading or trailing whitespace.
function result = parse_yes_no(str)
    str = strtrim(str);
    str = lower(str);
    
    switch str
        case 'yes'
            result = true;
        case 'true'
            result = true;
        case 'on'
            result = true;
        case 'no'
            result = false;
        case 'nope'
            result = false;
        case 'false'
            result = false;
        case 'off'
            result = false;
        otherwise
            error('Input not understood.')
    end
end

