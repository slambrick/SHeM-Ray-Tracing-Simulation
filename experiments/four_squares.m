% load scan data
dirname = '../results/four_squares_surf';
data = load(fullfile(dirname, 'scatteringData.mat'));

% extract counter values
data = data.simulationData;
counts = data.cntrSum{1};

% overall size of scan
xlims = [-0.5, 0.5];
zlims = [-0.5, 0.5];

% coords to slice one square
xslice = [0.051, 0.399];
zslice = [0.051, 0.399];

% make three slice of the background and compute counts
bck1 = slice_by_coords(counts, xlims, zlims, [-0.499, 0.499], [-0.499, -0.401]);
bck2 = slice_by_coords(counts, xlims, zlims, [-0.499, 0.499], [-0.049, +0.049]);
bck3 = slice_by_coords(counts, xlims, zlims, [-0.499, 0.499], [+0.401, +0.499]);
backgnd = [bck1; bck2; bck3];
[bckgnd_avg, bckgnd_err] = area_props(backgnd, 'Background');
square_counts_and_contrast(backgnd, bckgnd_avg, bckgnd_err, 'Background')

% slice each of the squares and compute contast to background
top_left = slice_by_coords(counts, xlims, zlims, -flip(xslice), -flip(zslice));
square_counts_and_contrast(top_left, bckgnd_avg, bckgnd_err, 'Top left');

top_right = slice_by_coords(counts, xlims, zlims, xslice, -flip(zslice));
square_counts_and_contrast(top_right, bckgnd_avg, bckgnd_err, 'Top right');

bot_right = slice_by_coords(counts, xlims, zlims, xslice, zslice);
square_counts_and_contrast(bot_right, bckgnd_avg, bckgnd_err, 'Bottom right');

bot_left = slice_by_coords(counts, xlims, zlims, -flip(xslice), zslice);
square_counts_and_contrast(bot_left, bckgnd_avg, bckgnd_err, 'Bottom left');

%% functions
function [avg, stddev] = area_props(counts, name)
    avg = mean(counts, 'all');
    stddev = std(counts, 0, 'all');
end

function [ctr, ctr_err] = mich_contrast(avg1, err1, avg2, err2)
    ctr = (avg1 - avg2) / (avg1 + avg2);
    ctr_err = 2 / (avg1 + avg2)^2 * sqrt((avg1 * err2)^2 + (avg2 * err1)^2);
end

function square_counts_and_contrast(counts, backgnd_avg, backgnd_err, name)
    [avg, stddev] = area_props(counts);
    fprintf('%s counts: %f +/- %f\n', name, avg, stddev);
    % disp([name ' counts: ' num2str(avg) ' +/- ' num2str(stddev)])

    [ctr, ctr_err] = mich_contrast(avg, stddev, backgnd_avg, backgnd_err);
    fprintf('%s contrast: %f +/- %f\n\n', name, ctr, ctr_err);
    % disp([name ' contrast: ' num2str(ctr) ' +/- ' num2str(ctr_err)])
end
