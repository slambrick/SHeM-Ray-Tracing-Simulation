CONST = 278.5085;

mass = [28, 197, 59, 195, 52];
% th_d = [625, 170, 375, 230, 460];
th_d = [230, 83, 220, 110, 175];
energy = 65;
inc_angle = 45;
temperature = 298;

w = CONST .* energy .* temperature ./ mass ./ (th_d .^ 2) * 2 ...
    .* (1 - cosd(180 - 2.*inc_angle));
dwf = exp(-w)

ctr_direct = tanh(CONST .* temperature .* energy .* (1 - cosd(180 - 2.*inc_angle))...
        .* (1/mass(1)/(th_d(1)^2) - 1./mass./(th_d.^2)))
ctr_from_dwf = (dwf - dwf(1)) ./ (dwf + dwf(1))
