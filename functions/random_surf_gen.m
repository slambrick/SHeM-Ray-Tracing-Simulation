function [f, xs] =  random_surf_gen(N, rL, h, cl)
% N - surface points
% rL - length of surface
% h - rms height
% cl - correlation length

h_RMS = h;
Dx = rL/(N - 1);
lambd = 0.5*cl^(2/3);

[f, xs] = random_surf_gen_core(h_RMS, Dx, lambd, 1, N);
end
        
