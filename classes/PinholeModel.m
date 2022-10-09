classdef PinholeModel
    % C class that defines the "N circle model" of detectors with a couple
    % of useful methods. The information needs to be transformed into a
    % struct to be passed into the C code.
    properties
        backwall_represent
        n_detectors
        circle_plate_radius
        aperture_axes
        aperture_c
        aperture_rotate
        in_use
        material
    end
    
    methods
        function obj = PinholeModel(material, represent, n_detectors, ...
                radius, axes, centres, aperture_rotate)
            if nargin == 1
                obj.backwall_represent = NaN;
                obj.n_detectors = NaN;
                obj.circle_plate_radius = NaN;
                obj.aperture_axes = NaN;
                obj.aperture_c = NaN;
                obj.aperture_rotate = NaN;
                obj.in_use = false;
                obj.material = material;
            elseif nargin == 7
                obj.backwall_represent = represent;
                obj.n_detectors = n_detectors;
                obj.circle_plate_radius = radius;
                obj.aperture_axes = axes;
                obj.aperture_c = centres;
                obj.aperture_rotate = aperture_rotate;
                obj.in_use = true;
                obj.material = material;
            else
                error('Only accepts 1 or 7 arguments');
            end
        end

        % Creates a plot of the detectors for this geometry
        function [f, ax] = plot_detectors(obj)
            f = figure;
            ax = axes(f);
            for i_=1:obj.n_detectors
                a = obj.aperture_axes(i_*2 - 1);
                b = obj.aperture_axes(i_*2);
                x0 = obj.aperture_c(i_*2 - 1);
                y0 = obj.aperture_c(i_*2);
                th = obj.aperture_rotate(i_);
                t=-pi:0.01:pi;
                x=a*cos(t);
                y=b*sin(t);
                tmp = x*cosd(th) - y*sind(th);
                y = x*sind(th) + y*cosd(th);
                x = tmp;
                x = x + x0;
                y = y + y0;
                plot(ax,x,y, 'DisplayName', ['Detector ' num2str(i_)]);
                hold on
            end
            plot(ax, [0], [0], 'o', 'MarkerFaceColor', 'red', 'Displayname', 'Design scattering point');
            xlabel(ax, 'x/mm')
            ylabel(ax, 'y/mm')
            legend(ax)
            hold off
            axis equal
            title('"N Circle" detector setup')
        end
        
        % NOTE: if the output names of this function are changed then changes
        % need to be made in the C code!
        function p = to_struct(obj)
            if obj.in_use
                p.plate_represent = obj.backwall_represent;
                p.n_detectors = obj.n_detectors;
                p.circle_plate_r = obj.circle_plate_radius;
                p.aperture_axes = obj.aperture_axes;
                p.aperture_c = obj.aperture_c;
                p.aperture_rotate = obj.aperture_rotate;
                p.material = obj.material;
            else
                error('The simple model of the pinhole plate has not been initialised');
            end
        end
    end
end

