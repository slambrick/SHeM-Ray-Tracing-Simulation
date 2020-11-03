close all
clear
clc

fpath = ['/home/sam/SMF_local/scatter_from_potential/non-useless-results/' ... 
    'roughness_investigation1_0060/'];
fname = [fpath, 'surface_used.csv'];

[xs, ys] = tracing2D.load_2Dsurface(fname);
xs = xs';
ys = ys';

disp('RMS height: ')
disp(std(ys))
[~, cl, ~] = roughSurf1D.acf(ys, xs);
disp('Correlation lenght: ')
disp(cl)

keyboard
n_rays = 20000;
init_angle = 20;

x_ray = linspace(-50, 50, n_rays);
y_ray = repmat(15, 1, n_rays);
init_pos = [x_ray; y_ray];
init_dir = [sind(init_angle); -cosd(init_angle)];

[thetas, num_scatters] = tracing2D.random_scatter('xs', xs, 'hs', ys, ...
    'n_rays', n_rays, 'scattering', 'specular', 'make_plots', false, 'init_pos', ...
    init_pos, 'init_dir', init_dir);

close all

figure
histogram(thetas, 'Normalization', 'pdf')
xlabel('\theta/^o')
ylabel('Probability')
xlim([-90, 90])
print([fpath, 'ray_tracing_hist.eps'], '-depsc')

save([fpath, 'ray_traced_results.mat'])