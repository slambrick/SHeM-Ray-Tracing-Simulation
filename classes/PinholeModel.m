classdef PinholeModel
    
    properties
        backwall_represent
        n_detectors
        circle_plate_radius
        aperture_axes
        aperture_c
        in_use
        material
    end
    
    methods
        function obj = PinholeModel(material, represent, n_detectors, radius, axes, centres)
            if nargin == 1
                obj.backwall_represent = NaN;
                obj.n_detectors = NaN;
                obj.circle_plate_radius = NaN;
                obj.aperture_axes = NaN;
                obj.aperture_c = NaN;
                obj.in_use = false;
                obj.material = material;
            elseif nargin == 6
                obj.backwall_represent = represent;
                obj.n_detectors = n_detectors;
                obj.circle_plate_radius = radius;
                obj.aperture_axes = axes;
                obj.aperture_c = centres;
                obj.in_use = true;
                obj.material = material;
            else
                error('Only accepts 1 or 6 arguments');
            end
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
                p.material = obj.material;
            else
                error('The simple model of the pinhole plate has not been initialised');
            end
        end
    end
end

