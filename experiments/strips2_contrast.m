sims_root = '../results';

% load data
sim_data = load(...
    fullfile(sims_root, 'dw_si_au_contrast', 'scatteringData.mat')).simulationData;

avg_a = zeros(1, 10);
avg_b = zeros(1, 10);

err_a = zeros(1, 10);
err_b = zeros(1, 10);

inc_energy = 20:8:92;

top_z = [0.501, 0.899];
bot_z = [0.101, 0.499];

for idx = 1:10
    x_slice = [0.011 + 0.1*(idx-1), 0.099*idx];

    strip_a = slice_by_coords(sim_data.cntrSum{1},...
        [0.0, 1.0], [0.0, 1.0], x_slice, bot_z);
    avg_a(idx) = mean(strip_a, 'all');
    err_a(idx) = std(strip_a, 0, 'all');

    strip_b = slice_by_coords(sim_data.cntrSum{1},...
        [0.0, 1.0], [0.0, 1.0], x_slice, top_z);
    avg_b(idx) = mean(strip_b, 'all');
    err_b(idx) = std(strip_b, 0, 'all');
end

contrast_val = abs(avg_a - avg_b) ./ (avg_a + avg_b);

contrast_err = sqrt((avg_a .* err_b) .^ 2 + (avg_b .* err_a) .^ 2) .* 2 ./ ((avg_a + avg_b) .^ 2);

figure;
errorbar(inc_energy, contrast_val, contrast_err, 'LineWidth', 1.5);
% xlim([0.05, 0.85])
xlabel('E_i / meV')
ylabel('Michelson Contrast')
set(gca, 'FontSize', 18)
set(gcf, 'Position', [100 100 800 600])

figure;
errorbar(inc_energy, log(avg_a), err_a./avg_a, 'LineWidth', 1.5);
hold on

errorbar(inc_energy, log(avg_b), err_b./avg_b, 'LineWidth', 1.5);

legend('Gold', 'Silicon')

xlabel('E_i / meV')
ylabel('log(I / a.u.)')
set(gca, 'FontSize', 18)
set(gcf, 'Position', [100 100 800 600])
