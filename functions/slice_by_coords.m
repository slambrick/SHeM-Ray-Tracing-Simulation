function B = slice_by_coords(A, xlims, ylims, xslice, yslice)
    %   Function that will slice an image (i.e. a 2D array A)
    %   that spans physical coordinates xlims, ylims
    %   to the part contained between coordinates xslice, yslice
    %   INPUTS:
    %       - A = the input 2D array
    %       - xlims, ylims = arrays of 2 doubles each giving
    %         the coordinate limits of the image
    %       - xslice, yslice = arrays of 2 doubles each
    %         giving the coordiante limits of the desired slice

    dims = size(A);
    nrows = dims(1); ncols = dims(2);
    if length(dims) ~= 2
        error("A must be 2-dimensional array")
    end

    % the physical size of the sample
    x_size = xlims(2) - xlims(1);
    y_size = ylims(2) - ylims(1);

    % the percentage along the x axis where the slice starts and ends
    x_p1 = (xslice(1) - xlims(1)) / x_size;
    x_p2 = (xslice(2) - xlims(1)) / x_size;

    % percentage along the y axis where slice starts and ends
    y_p1 = (yslice(1) - ylims(1)) / y_size;
    y_p2 = (yslice(2) - ylims(1)) / y_size;

    % turn percentages into array indices
    col1 = ceil(x_p1*(ncols-1)) + 1;
    col2 = floor(x_p2*(ncols-1)) + 1;

    row1 = ceil(y_p1*(nrows-1)) + 1;
    row2 = floor(y_p2*(nrows-1)) + 1;

    B = A(row1:row2, col1:col2);
