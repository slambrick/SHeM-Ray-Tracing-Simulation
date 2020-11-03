classdef BackWall
    %BACKWALL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        represent_wall
        n_detectors
        circle_plate_r
        aperture_axes
        aperture_c
        material
    end
    
    methods
        function obj = Sphere(make, c, r, material)
            if nargin == 1
                c = [0,0,0];
                r = 0;
                material = 'none';
            end
            obj.centre = c;
            obj.radius = r;
            obj.make = make;
            obj.material = material;
        end
        
        % NOTE: if the output names of this function are changed then changes
        % need to be made in the C code!
        function s = to_struct(obj)
            s.c = obj.centre;
            s.r = obj.radius;
            s.make = obj.make;
            s.material = obj.material;
        end
    end
end
