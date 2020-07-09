% bulk temps
% obs_contrast = [0, -0.387, -0.145, -0.028, 0.003];
% sig_contrast = [0, 0.009, 0.008, 0.008, 0.008];

% surface temps
obs_contrast = [0, 0.0, 0.27, 0.14, 0.0];
sig_contrast = [0, 0.03, 0.03, 0.03, 0.03];

CONST = 278.5085;

mass = [28, 197, 59, 195, 52];
% th_d = [625, 170, 375, 230, 460];   % bulk
th_d = [230, 83, 220, 110, 175];  % surface
energy = 65;
inc_angle = 45;
temperature = 298;

w = CONST .* energy .* temperature ./ mass ./ (th_d .^ 2) * 2 ...
    .* (1 - cosd(180 - 2.*inc_angle));
dwf = exp(-w);

dwf_sum = dwf + dwf(1);
dwf_dif = dwf - dwf(1);

diffuse_contrib = (dwf_dif - obs_contrast .* dwf_sum) ./...
    (obs_contrast .* (2 - dwf_sum) + dwf_dif);

diffuse_contrib_err = abs(2 .* sig_contrast .* dwf_dif ./...
    ((obs_contrast .* (2 - dwf_sum) + dwf_dif).^2));

% using the values for the diffuse contrib obtained above,
% can calculate mean and error diffuse contribution
values = [0.0285, 0.0281, 0.0253, 0.0233, 0.0219];
errs   = [0.0043, 0.0070, 0.0064, 0.0182, 0.1051];
weigh  = 1 ./ (errs .^ 2);

format long
average_diffuse_contr = sum(values .* weigh) / sum(weigh)
error_diffuse_contr   = sqrt(sum(weigh .* (values - average_diffuse_contr).^2)...
    ./ sum(weigh) / (length(values) - 2))
