% Copyright (c) 2018-20, Sam Lambrick.
% All rights reserved.
% This file is part of the Sub-beam Ray Tracing Simulation, subject to the  
% GNU/GPL-3.0-or-later.
%
% A script for plotting the weighted histograms of the results from
% 'loop_roughness.m', plots the weighted histogram of the final angles along
% with a crude single parameter fit to the broad specular distribution (+ cosine
% distribution and the specular location).
clear
close all
clc

%% Parameters

% The data file produced by 'loop_roughnes.m'
data_file = '../results1D/increasing_roughness011.mat';

% Plot example profiles of the surface scattered off
profiles = false;

% Directory to save the plots to
the_path = '../results1D/';

% Incidence angle
init_angle = 45;

%% Preperation

loadpath

load(data_file)

% Plot of the proportion of multiple scattering
[f, proportion_multiple] = tracing2D.multiple_scattering_hist(num_scatters, ratios);
print([the_path 'proportion_multiple_scattering.eps'], '-depsc')

% Histogram plots
generate_histograms(ratios, thetas, num_scatters, 'both', the_path, init_angle);
generate_histograms(ratios, thetas, num_scatters, 'single', [the_path 'single/'], init_angle);
generate_histograms(ratios, thetas, num_scatters, 'multiple', [the_path 'multiple/'], init_angle);

%% Gnerate histograms for one type of scattering
function generate_histograms(ratios, thetas, num_scatters, scattering, the_path, init_angle)

    if ~exist(the_path, 'dir')
        mkdir(the_path)
    end
    if ~exist([the_path 'pngs'], 'dir')
        mkdir([the_path 'pngs']);
    end

    % The cosine distribution, this will be added to all plots
    normalisation = integral(@(x) cosd(x), -90, 90);
    xs_cos = -90:0.1:90;
    ys_cos = cosd(xs_cos)/normalisation;

    % The relation between sigma and sigmaPs, we crudly estimate the appropriate
    % broad specular by numerically aquiring the parameter sigma from the standard
    % devaition of the histogram.
    sigma_Ps = 0.1:0.1:360;
    sigma_Ps = sigma_Ps*pi/180;
    sigma_s = zeros(size(sigma_Ps));
    for i_=1:length(sigma_Ps)
        [~, sigma_s(i_)] = sigma_of_sigmaP(sigma_Ps(i_), init_angle*pi/180);
    end
    sigma_s = sigma_s*180/pi;
    sigma_Ps = sigma_Ps*180/pi;

    for i_=1:length(ratios)
        %% Sort by scattering
        switch scattering
            case 'single'
                ind = num_scatters(i_,:) == 1;
                thetas_use = thetas(i_,:);
                thetas_use = thetas_use(ind);
                profiles = false;
            case 'multiple'
                ind = num_scatters(i_,:) > 1;
                thetas_use = thetas(i_,:);
                thetas_use = thetas_use(ind);
                profiles = false;
            case 'both'
                thetas_use = thetas(i_,:);
                profiles = true;
            otherwise
                error('Specify single, multiple, or both for scattering')
        end

        %% Plot of the resultant scattering distribution
        if profiles
            figure('units', 'centimeters', 'Position', [0 0 14*1.5 15*1.5]);
            subplot(3, 1, 1:2)
        else
            figure('units', 'centimeters', 'Position', [0 0 14*1.5 10*1.5]);
        end
        the_hist = histogram(thetas_use, 30, 'Normalization', 'pdf', ...
            'LineStyle', 'none', 'FaceAlpha', 0.7, 'FaceColor', [0.2 0.1 0.7]);
        title(['RMS height/correlation length = ' num2str(ratios(i_))]);
        hold on

        % Plot a smoothed version of the data as a pdf
        if ~isempty(thetas_use)
            pd = fitdist(thetas_use', 'Kernel', 'Bandwidth', the_hist.BinWidth);
            y = pdf(pd, -90:0.1:90);
            plot(-90:0.1:90, y, 'Color', [0 0 1], 'LineWidth', 1);
        end

        % Plot a cosine
        plot(xs_cos, ys_cos, 'r')

        % Plot my broad specular distribution
        % First find the parameter standard deviation that corresponds to the
        % standard deviation of the 
        the_sigma = std(the_hist.Data);
        if ~isempty(thetas_use)
            if the_sigma >= max(sigma_s)
                % The broadest the model can represent
                parameter_sigma = 180;
            else
                ind = sigma_s < the_sigma;
                ind = sum(ind);
                parameter_sigma = sigma_Ps(ind);
            end

            if ratios(i_) <= 0.32
                % Forward scattering
                ys_spec = broad_specular_distribution(xs_cos, -init_angle, parameter_sigma);
            else
                % Back scattering
                ys_spec = broad_specular_distribution(xs_cos, init_angle, parameter_sigma);
            end
            plot(xs_cos, ys_spec, 'g')
        else
            the_sigma = NaN;
            parameter_sigma = NaN;
            ys_spec = 0.012;
        end
        fprintf('sigma = %f\n', the_sigma)
        fprintf('sigma_P = %f\n\n', parameter_sigma)

        % Emphasis where the specular is
        plot([init_angle init_angle], [0 1.2*max(ys_spec)], 'k', 'Linewidth', 1.5)

        % Add a legend
        if length(thetas_use > 0)
            if ratios(i_) <= 0.32
                legend('Histogram', 'Smoothed pdf', 'Cosine', 'Broad specular', ...
                    'Specular', 'Location', 'NorthWest')
            else
                legend('Histogram', 'Smoothed pdf', 'Cosine', 'Backscattering', ...
                    'Specular', 'Location', 'NorthEast')
            end
        else
            legend('Histogram', 'Cosine', 'Specular', 'Location', 'NorthWest')
        end
        hold off
        xlabel('\theta/\circ')
        ylabel('P(\theta)')
        xlim([-90, 90])
        ylim([0, 1.2*max(ys_spec)]);

        %% Plot of a sample of the random surface
        if profiles
            subplot(3, 1, 3)
            [hs, xs] = rsgene1D(5000, 50, ratios(i_), 1);
            plot(xs, hs)
            axis equal
            xlabel('x')
            ylabel('Height')
            title(['RMS height/correlation length = ' num2str(ratios(i_))])
        end

        if mod(ratios(i_), 0.1) == 0
            print(gcf, [the_path 'hist_ratio_' num2str(ratios(i_)) '0.eps'], '-depsc')
            print(gcf, [the_path 'pngs/hist_ratio_' num2str(ratios(i_)) '0.png'], '-dpng', '-r0')
        else
            print(gcf, [the_path 'hist_ratio_' num2str(ratios(i_)) '.eps'], '-depsc')
            print(gcf, [the_path 'pngs/hist_ratio_' num2str(ratios(i_)) '.png'], '-dpng', '-r0')
        end
    end
end