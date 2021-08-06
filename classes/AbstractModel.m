classdef AbstractModel
    
    properties
        theta
        phi
        direction
        half_cone_angle
        in_use
    end
    
    methods
        function obj = AbstractModel(theta, phi, half_cone_angle)
            if nargin == 0
                obj.theta = NaN;
                obj.phi = NaN;
                obj.direction = NaN;
                obj.half_cone_angle = NaN;
                obj.in_use = false;
            else
                obj.theta = theta;
                obj.phi = phi;
                obj.direction = [cosd(phi)*sind(theta), cosd(theta), sind(phi)*sind(theta)];
                obj.half_cone_angle = half_cone_angle;
                obj.in_use = true;
            end
        end
        
        function p = to_struct(obj)
            if obj.in_use
                p.direction = obj.direction;
                p.half_cone_angle = obj.half_cone_angle*pi/180;
            else
                error('The simple model of the pinhole plate has not been initialised');
            end
        end
    end
end

