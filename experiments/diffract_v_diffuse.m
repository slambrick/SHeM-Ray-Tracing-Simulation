dir1 = '../results/peaks_scan';
dir2 = '../results/peaks_diffuse';
dirname = '../results/diffract_v_diffuse';
if ~exist(dirname, 'dir')
    mkdir(dirname);
end

simdata1 = load(fullfile(dir1, 'scatteringData.mat')).simulationData{1};
simdata2 = load(fullfile(dir2, 'scatteringData.mat')).simulationData;

data_arr = {simdata1, simdata2};

n_rays = data_arr{1}.rays_per_pixel;
min_count = n_rays*2;
max_count = -1;

% loop through scans to find min and max brightness values to scale
for idx = 1:2
    counters = data_arr{idx}.cntrSum{1};

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
fnames = {'diffract.png', 'diffuse.png'};

for idx = 1:2
    sim_data = data_arr{idx};

    fname = fnames{idx};
    fname = fullfile(dirname, fname);

    img = sim_data.imageAll('scale', 'manual', 'specifyScale', limits);
    set(gca, 'xtick', [], 'ytick', []);
    xlabel(''); ylabel('');
    exportgraphics(gca, fname, 'Resolution', 300, 'BackgroundColor', 'black');
end
