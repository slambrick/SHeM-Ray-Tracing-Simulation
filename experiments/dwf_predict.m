% bulk temps
% obs_contrast = [0, -0.387, -0.145, -0.028, 0.003];
% sig_contrast = [0, 0.009, 0.008, 0.008, 0.008];

% the diffuse contribution
f = 0.0275;
df = 0.0009;

% surface temps
obs_contrast = [0.015, 0.005, 0.001, -0.009];
sig_contrast = [0.005, 0.003, 0.003,  0.004];

CONST = 278.5085;

mass = [28, 197, 59, 195, 52];
% th_d = [625, 170, 375, 230, 460];   % bulk
th_d = [230, 83, 220, 110, 175];  % surface
energy = 63;
inc_angle = 45;
temperature = 294;

w = CONST .* energy .* temperature ./ mass ./ (th_d .^ 2) * 2 ...
    .* (1 - cosd(180 - 2.*inc_angle));
dwf = exp(-w);

% assuming the first one is correct.
ref_dwf = dwf(1);
% get the DWF from contrast expression
dwf_predicted = ref_dwf .* (1 + obs_contrast) ./ (1 - obs_contrast) ...
    + 2 * f / (1-f) .* obs_contrast ./ (1 - obs_contrast);
log_dwf = log(dwf_predicted);

% then the debye temp th_d from the DWF prediction
th_d_predicted = sqrt(...
    - 4 * CONST * temperature * energy * (cos(inc_angle))^2 ...
    ./ (mass(2:end) .* log_dwf))


% for the errors, calculate the contributions first:
ctr_error = sig_contrast .* (ref_dwf * 2 ./ ((1 - obs_contrast).^2)...
    + 2 * f * (1-f) ./ (((1-f) .* (1-obs_contrast)).^2));
% f_error = df * 2 .* obs_contrast .* (1+obs_contrast) ./ (((1-f) * (1+obs_contrast)).^2);
f_error = 0;
dwf_error = sqrt(ctr_error.^2 + f_error.^2);
th_d_error = dwf_error ./ dwf_predicted .* th_d_predicted


%% Check
contr_check = (1-f)*(dwf_predicted - ref_dwf) ./...
    (2*f + (1-f) * (dwf_predicted + ref_dwf))
