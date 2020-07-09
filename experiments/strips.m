sims_root = '../results';

% load data
sim_data = load(...
    fullfile(sims_root, 'strips/scatteringData.mat')).simulationData;

averages = zeros(1, 8);
errors = zeros(1, 8);
spreads = 0.1:0.1:0.8;

for idx = 1:1:8
    x_slice = [0.1*idx + 0.005, 0.1*(idx+1) - 0.005];
    strip = slice_by_coords(sim_data.cntrSum{1},...
        [0.0, 1.0], [0.0, 1.0], x_slice, [0.101, 0.899]);
    averages(idx) = mean(strip, 'all');
    errors(idx) = std(strip, 0, 'all');
end

bck = slice_by_coords(sim_data.cntrSum{1},...
    [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 0.099]);
bck_level = mean(bck, 'all');

x_axis_lims = [0.05, 0.9];

y_axis_scale = 1e3;

errorbar(spreads, averages/y_axis_scale, errors/y_axis_scale, 'LineWidth', 2, 'CapSize', 12);
hold on
plot(x_axis_lims, [bck_level, bck_level]/y_axis_scale, 'k--', 'LineWidth', 2);

xlim(x_axis_lims);
xlabel('\sigma (rad)');
ylabel('Average Counts / 10^{3}');
set(gca, 'FontSize', 30);
set(gcf, 'Position', [100 100 700 600]);
