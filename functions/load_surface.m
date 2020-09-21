% Loads a previously generated surface for a text file for use in 2D ray tracing
function [x, y, n] = load_surface(fname)
    opts = delimitedTextImportOptions("NumVariables", 2);

    % Specify range and delimiter
    opts.DataLines = [10, Inf];
    opts.Delimiter = ",";

    % Specify column names and types
    opts.VariableNames = ["x", "y"];
    opts.VariableTypes = ["double", "double"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    % Import the data
    surfaceused = readtable(fname, opts);
    x = surfaceused.x;
    y = surfaceused.y;
    dx = x(2) - x(1);
    n = diff(y)/dx;
end

