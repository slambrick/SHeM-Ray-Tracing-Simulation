% Load the LiF movie dataset
dirname = '../results/LiF_movie';
data = load(fullfile(dirname, 'scatteringData.mat'));
dirname = fullfile(dirname, 'movie1');
if ~exist(dirname, 'dir')
    mkdir(dirname);
end

% what are the y displacements
range_y = data.scan_inputs.range1D;
raster_movement_y = data.scan_inputs.raster_movment1D;
ys = range_y(1):raster_movement_y:range_y(2);

% extract the scan data
simulationData = data.simulationData;
dist_to_sample = simulationData{1}.dist_to_sample;

% actual y distances = displacement + original dist to sample
ys = ys + dist_to_sample;

% check that all points were recorded
if length(simulationData) ~= length(ys)
    error(['There must be as many y points as scans. Currently ' length(simulationData) ...
            ' scans and ' length(ys) ' y points.']);
end

% how many rays did we use
n_rays = simulationData{1}.rays_per_pixel;
min_count = n_rays*2;
max_count = -1;

% loop through scans to find min and max brightness values to scale
for idx = 1:length(ys)
    counters = simulationData{idx}.cntrSum{1};

    current_max = max(counters, [], 'all');
    current_min = min(counters, [], 'all');

    if current_max > max_count
        max_count = current_max;
    end
    if current_min < min_count
        min_count = current_min;
    end
end

limits = [min_count, max_count]

% loop through scans and produce figures
for idx = 1:length(ys)
    y = ys(idx);
    sim_data = simulationData{idx};

    fname = sprintf('z%.2fmm.png', y);
    fname = fullfile(dirname, fname);

    img = sim_data.imageAll('scale', 'manual', 'specifyScale', limits);
    set(gca, 'xtick', [], 'ytick', []);
    xlabel(''); ylabel('');
    text(-0.15, -0.13, sprintf('z = %.2f mm', y), 'Color', 'white', 'FontSize', 7);
    % imwrite(img, fname);
    exportgraphics(gca, fname, 'Resolution', 300, 'BackgroundColor', 'none');
    % close(gcf);
end

% close all figures
close all
