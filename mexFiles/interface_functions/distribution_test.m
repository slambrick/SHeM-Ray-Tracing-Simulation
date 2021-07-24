% Copyright (c) 2020-21, Dan Seremet, Sam Lambrick.
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
    % The lattice needs to be set if the diffraction models are to be used
    % the default is for a square lattice
    lattice = [0.0 1.0 0.0 0.0 0.0 1.0];
    
    switch scattering
        case 'diffraction'
            % A diffraction model that provides the lattice vectors
            diff_params = [
                6,   6,   0.1996,... % maxp and maxq, lambda/a
                cosd(14), sind(14), -sind(14), cosd(14),... % rec lattice vectors
                0.0316,  2];       % peak sigma and envelope sigma

            material.function = 'diffraction';
            material.params = [0.6, diff_params];
            material.color = [0.8, 0.8, 1.0];
        case 'cosine'
            % A purely diffuse distribution following a cosine function
            % centred on the surface normal.
            material.function = 'cosine';
            material.params = [];
            material.color = [1.0, 0.8, 0.8];
        case 'uniform'
            % A completely uniform distribution with equal probability of
            % scattering in all directions
            material.function = 'uniform';
            material.params = [45]; % Need a cutoff to stop scattering into 
                                    % the sample plane, can't allow 90deg
            material.color = [1.0, 0.3, 0.3];
        case 'broad specular'
            % A broadened specular distibution based on a Gaussian, a
            % proportion of the scattering can be diffuse
            material.function = 'broad_specular';
            material.params = [0.5, 50];   % 50% specular with S.D. of 50deg
            material.color = [0.3, 1, 0.3];
        case 'cosine specular'
            % A cosine scattering model centred on the specular condition.
            material.function = 'cosine_specular';
            material.params = [];
            material.color = [0.8, 1, 0.8];
        case 'diffraction2'
            % A diffraction model that uses the lattice vectors described
            % in 3D in the sample object
            material.function = 'diffraction2';
            material.params = [0.6, ...          % Proportion that is diffuse
                               6, 6, 0.1996, ... % maxp and maxq, lambda/a
                               0.0316, 2];       % peak sigma and envelope sigma
            material.colour = [0.5, 0.5, 1.0];
        otherwise
            error('Specified type of scattering not recognised');
    end

    [theta, phi, final_dir] = distributionTestMex(num_rays, direction, material, normal, lattice);

    %[theta, phi] = calc_angles(final_dir', normal', direction');
    
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

