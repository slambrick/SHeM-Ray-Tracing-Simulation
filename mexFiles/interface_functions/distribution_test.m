% Copyright (c) 2020, Dan Seremet, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the
% GNU/GPL-3.0-or-later.
%
% Interface for testing the generation of scattering distributions
function [theta, phi] = distribution_test(varargin)
    
    for i_=1:2:length(varargin)
        switch varargin{i_}
            case 'n_rays'
                num_rays = varargin{i_+1};
            case 'init_dir'
                direction = varargin{i_+1};
            case 'scattering'
                scattering = varargin{i_+1};
            case 'recompile'
                recompile = varargin{i_+1};
            otherwise
                warning(['Input' num2str(i_) 'not recognised']);
        end
    end
    
    if ~exist('recompile', 'var')
        recompile = false;
    end
    
    mexCompile(recompile);
    
    if ~exist('num_rays', 'var')
        num_rays = 1e6;
    end
    
    %direction = [sin(pi/6), -cos(pi/6), 0];
    direction = direction/norm(direction);
    normal = [0, 1, 0];
    
    switch scattering
        case 'diffraction'
            diff_params = [
                6,   6,   0.1996,... % maxp and maxq, lambda/a
                cosd(14), sind(14), -sind(14), cosd(14),... % rec lattice vectors
                0.0316,  2];       % peak sigma and envelope sigma

            material.function = 'diffraction';
            material.params = [0.6, diff_params];
            material.color = [0.8, 0.8, 1.0];
        case 'cosine'
            material.function = 'cosine';
            material.params = [];
            material.color = [1.0, 0.8, 0.8];
        case 'uniform'
            material.function = 'uniform';
            material.params = [];
            material.color = [1.0, 0.3, 0.3];
        case 'broad specular'
            material.function = 'broad_specular';
            % 50% specular with S.D. of 20deg
            material.params = [0.5, 20*pi/180];
            material.color = [0.3, 1, 0.3];
        case 'cosine specular'
            material.function = 'cosine_Specular';
            material.params = [];
            material.color = [0.8, 1, 0.8];
        otherwise
            error('Specified type of scattering not recognised');
    end

    [theta, phi] = distributionTestMex(num_rays, direction, material, normal);

    plot_distribution_3d(sin(theta), phi, 1, 100, '\theta')
    plot_distribution_slice(theta, phi, 0, 0.09, 200)
    
    
    % Polar density plot of the generated values
    figure
    [n, c] = hist3([theta*180/pi, phi*180/pi + 180], {2.5:5:87.5, 2.5:5:357.5}, 'Normalization', 'pdf');
    n2 = zeros(size(n));
    for i_=1:length(c{1})
        if c{1}(i_) ~= 0
            n2(i_,:) = n(i_,:)/sind(c{1}(i_));
        else
            n2(i_,:) = 0;
        end
    end
    n2 = n2/max(n2(:));
    polarPcolor(90*(c{1} - 2.5)/max(c{1} - 2.5), 360*(c{2} - 2.5)/max(c{2} - 2.5), n2', 'colBar', false)
    ylabel('Relative intensity')
    print(gcf, 'density_plot_generated60.eps', '-depsc')
    title('Generated directions')
end
