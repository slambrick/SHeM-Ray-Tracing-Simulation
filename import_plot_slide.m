function [I, zs, det_pos] = import_plot_slide(dfile)
load(dfile);

m_u = 1.6605e-27;
m_He = 4*m_u; 
k_B = 1.380649e-23;
h = 6.62607015e-34;
E = 64e-3*1.6e-19;
k = 2*pi*sqrt(2*m_He*E)/h;

n_slide = length(simulationData);
n_z = length(simulationData{1}.sample_positions);
I = zeros(n_slide, n_z);
zs = simulationData{1}.sample_positions;
det_pos = -3.5:0.1:3.5;
for i_=1:n_slide
    s = simulationData{i_};
    I(i_,:) = s.single_scattering + s.multiple_scattering;
end

zs2 = repmat(zs, [length(det_pos), 1]);
det_pos2 = repmat(det_pos', [1, length(zs)]);

a = 7;
DKx = k*(sin(atan((a - zs2)./zs2)) - 1/sqrt(2));
DKy = k*sin(atan(det_pos2./zs2));

figure
s = pcolor(DKx*1e-9, DKy*1e-9, I);
s.LineWidth = 0;
xlabel('\Delta K_x/nm^{-1}')
ylabel('\Delta K_y/nm^{-1}')
axis equal
end