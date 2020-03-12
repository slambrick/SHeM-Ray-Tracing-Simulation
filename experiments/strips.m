sims_root = '../results';

% load data
sim_data = load(...
    fullfile(sims_root, 'strips/scatteringData.mat')).simulationData;

averages = zeros(1, 8);
errors = zeros(1, 8);
spreads = 0.1:0.1:0.8;

for idx = 1:1:8
    x_slice = [0.1*idx, 0.1*(idx+1)];
    strip = slice_by_coords(sim_data.cntrSum{1},...
        [0.0, 1.0], [0.0, 1.0], x_slice, [0.1, 0.9]);
    averages(idx) = mean(strip, 'all');
    errors(idx) = std(strip, 0, 'all');
end

errorbar(spreads, averages, errors, 'LineWidth', 1.5);
xlim([0.05, 0.85])
xlabel('\Delta\theta (rad)');
ylabel('Average Counts')
set(gca, 'FontSize', 18);
set(gcf, 'Position', [100 100 500 400])
