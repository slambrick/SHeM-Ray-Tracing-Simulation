% classdef BackWall
%     % May be depreciated, I need to check!
%     
%     properties
%         represent_wall
%         n_detectors
%         circle_plate_r
%         aperture_axes
%         aperture_c
%         material
%     end
%     
%     methods
%         function obj = Sphere(make, c, r, material)
%             if nargin == 1
%                 c = [0,0,0];
%                 r = 0;
%                 material = 'none';
%             end
%             obj.centre = c;
%             obj.radius = r;
%             obj.make = make;
%             obj.material = material;
%         end
% 
%         % Creates a plot of the detectors for this geometry
%         function [f, ax] = plot_detectors(obj)
%             f = figure;
%             ax = axes(f);
%             for i_=1:obj.n_detectors
%                 a = obj.aperture_axes(i_*2 - 1);
%                 b = obj.aperture_axes(i_*2);
%                 x0 = obj.aperture_c(i_*2 - 1);
%                 y0 = obj.aperture_c(i_*2);
%                 t=-pi:0.01:pi;
%                 x=x0+a*cos(t);
%                 y=y0+b*sin(t);
%                 plot(ax,x,y)
%             end
%         end
%         
%         % NOTE: if the output names of this function are changed then changes
%         % need to be made in the C code!
%         function s = to_struct(obj)
%             s.c = obj.centre;
%             s.r = obj.radius;
%             s.make = obj.make;
%             s.material = obj.material;
%         end
%     end
% end
