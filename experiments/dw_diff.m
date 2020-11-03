compile_tests()

%% General input parameters
num_rays = 1e6;
direction = [sin(pi/6), -cos(pi/6), 0];
direction = direction/norm(direction);
normal = [0, 1, 0];

diff_params = [
        6,   6,   0.1996,... % maxp and maxq
        1,   0,   0,   1,... % rec lattice vectors
        0.0316,  100];       % peak sigma and envelope sigma

dw_params = [20, 197, 298, 170, 0];

mat_dw.function = 'dw_specular';
mat_dw.params = [dw_params, 0.1];
mat_dw.color = [0.8, 0.8, 1.0];

[theta_1, phi_1] = distribution_test_mex(num_rays, direction, mat_dw, normal);

mat_spec.function = 'broad_specular';
mat_spec.params = [0, 0.1];
mat_spec.color = [0.8, 0.8, 1.0];

[theta_2, phi_2] = distribution_test_mex(num_rays, direction, mat_spec, normal);

%% set up histogramming
rad_limit = 1;
nbins = 100;

bin_values_1 = norm_histogram2(sin(theta_1), phi_1, rad_limit, nbins);
bin_values_2 = norm_histogram2(sin(theta_2), phi_2, rad_limit, nbins);
close(gcf);

%% plot the two original distributions
% plot_distribution_3d(sin(theta_1), phi_1, 1, 100, '\theta');
% plot_distribution_3d(sin(theta_2), phi_2, 1, 100, '\theta');

%% plot the differences
% first diff the two histograms
bin_values_diff = bin_values_1 - bin_values_2;

% find the x and y bin positions
bin_step = 2*rad_limit / nbins;
bin_centres = -(rad_limit - bin_step/2) : bin_step : (rad_limit - bin_step/2);

figure;
surf(bin_centres, bin_centres, bin_values_diff', 'EdgeColor', 'None');

colormap(custom_colormap(bin_values_diff, 0));

%% plot settings
view([0,0,1]);
colorbar;
xlim([-rad_limit, rad_limit])
ylim([-rad_limit, rad_limit])
% xlabel([name ' cos(\phi)'])
% ylabel([name ' sin(\phi)'])
yticks(-1:0.5:1);

xlabel('n_{f, x}')
ylabel('n_{f, y}')
zlabel('Difference in counts')
% title([name ' spatial distribution'])
set(gca, 'FontSize', 30)
set(gcf, 'Position', [100 100 900 800])
daspect([1 1 1])

function bin_values = norm_histogram2(rad, phi, rad_limit, nbins)
    rad_step = 2*rad_limit / nbins;
    rad_edges = -rad_limit:rad_step:rad_limit;

    pos2d = [rad.*cos(phi), rad.*sin(phi)];
    hist_1 = histogram2(pos2d(:, 1), pos2d(:, 2), rad_edges, rad_edges, 'FaceColor', 'flat');
    bin_values = hist_1.Values;
    % bin_values = hist_1.Values ./ max(hist_1.Values, [], 'all');
end