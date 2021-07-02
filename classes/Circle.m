classdef Circle
    
    properties
        centre
        radius
        normal
        make
        material
    end
    
    methods
        function obj = Circle(make, material, c, r, n)
            if nargin == 2
                c = [0,0,0];
                r = 0;
                n = [0,1,0];
            end
            obj.centre = c;
            obj.radius = r;
            obj.normal = n;
            obj.make = make;
            obj.material = material;
        end
        
        % NOTE: if the output names of this function are changed then changes
        % need to be made in the C code!
        function c = to_struct(obj)
            c.c = obj.centre;
            c.r = obj.radius;
            c.n = obj.normal;
            c.make = obj.make;
            c.material = obj.material;
        end
        
        function obj = move(obj, vect)
            obj.centre = obj.centre + vect;
        end
        
        function obj = scale(obj, fct)
            obj.radius = obj.radius*fct;
        end
    end
end

