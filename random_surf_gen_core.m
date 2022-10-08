function [f, xs] = random_surf_gen_core(h, Dx, lambd, s, N)
% Core routine to generate a random surface. Arguments should be
% chosen carefully. The number of points must be odd.

    if nargin == 4
        N = 10001;
    end
    if mod(N, 2) == 0
        error("Must be an odd number of points in the surface");
    end
    Z = normrnd(0.0, s, 1, N);
    ms = linspace(-round(N/2), round(N/2), N);
    e = exp(-abs(ms)*Dx/lambd);
    f = h*conv(Z, e, 'same');
    xs = ms*Dx;
end