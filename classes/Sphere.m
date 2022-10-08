classdef Sphere
    
    properties
        centre
        radius
        make
        material
        n_sphere
    end
    
    methods
        function obj = Sphere(make, material, c, r)
            if nargin == 2
                c = [0;0;0];
                r = 0;
            end
            obj.n_sphere = length(r);
            if length(r) ~= size(c, 2)
                error('You made a bad list of spheres');
            end
            obj.centre = c;
            obj.radius = r;
            obj.make = make;
            obj.material = material;
        end
        
        % NOTE: if the output names of this function are changed then changes
        % need to be made in the C code!
        function [s, n_sphere] = to_struct(obj)
            s.c = obj.centre(:)';
            s.r = obj.radius;
            s.make = obj.make;
            s.material = obj.material;
            n_sphere = obj.n_sphere;
        end
        
        function obj = move(obj, vect)
            obj.centre = obj.centre + vect;
        end
        
        function obj = scale(obj, fct)
            obj.radius = obj.radius*fct;
        end
    end
end

