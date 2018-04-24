classdef TriagSurface < handle
% Copyright (c) 2018, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
%
% TriagSurface.m 
%
% Contains a triangulated mesh surface, consisting of a series of vertices
% and triangles. Each triangle consists of three vertices and associated
% with each triangle is a normal direction a marker that determines the
% nature of the scattering off that region of sample.
%
% PROPERTIES:
%  vertices    - A 3xn list of coordinates for each vertex in the triangulation.
%  faces       - A 3xm list of indices giving which vertices make up triangles
%                in the surface mesh.
%  normals     - A 3xm list of the unit normals to each of the triangles
%                specifed by the faces property.
%  composition - A 1xm list specifying what form of scattering occurs off of
%                each triangle. Can be a value between 0 and 1 or can be 2.
%  nTriag      - The number of triangles in the surface mesh.
%  nVertices   - The number of veritces in the surface mesh.
    
    properties (SetAccess = private)
        vertices    % Vertices in the triangulation
        faces       % Contains the reference the vertices that make up the 
                    % triangles
	    normals     % Normals to the triangles
        composition % The amount of diffuse scattering for each triangle
        nTriag      % The number of triangles in the surface
        nVertices   % The number of vertices in the surface
    end % End properties
    
    methods
        function obj = TriagSurface(V, N, F, C)
        % Constructor for TriagSurface class. Can be called with 0 or 4
        % arguments. 0 arguments creats an empty object.
        %
        % INPUTS:
        %  V - The coordinates of the vertices.
        %  N - The normals to the faces of the triangles.
        %  F - The 'faces', contains the index of the vertices that make up
        %      each face.
        %  C - A marker determining the kind of scattering that occurs on
        %      each triangle.
        %
        % OUTPUT:
        %  obj - A TriagSurface object that contains the surface.
            if (nargin == 0)
                % Default, creates an empty object
            elseif (nargin == 4)
                if (size(size(V), 2) ~= 2)
                    error('Input arguments must all be 2D arrays')
                end
                if (size(V,2) ~= 3 && size(N,2) ~= 3 && size(F,2) ~= 2)
                    error('V, N, F inputs must all have three columns')
                end
                if (size(C,2) ~= 1)
                    error('C input must be a 1D column vector')
                end
                if (size(V,1) < 3)
                    error('There must be at least three vertices in the surface')
                end
                if (size(N, 1) ~= size(F,1) && size(N,1) ~= length(C))
                    error('Inputs F, N, C must have the same length')
                end
                if (sum(~(C <= 2)) ~= 0 || sum(~(C >= 0)))
                    error('The diffusivity level can only be betwen 0 and 2')
                end
                obj.vertices = V;
                obj.normals = N;
                obj.faces = F;
                obj.composition = C;
                obj.nTriag = length(C);
                obj.nVertices = size(V, 1);
            else
                error('Wrong number of input arguments');
            end
        end % End constructor
        
        function newobj = copy(obj)
        % Copys the object and returns the copy.
            newobj = TriagSurface(obj.vertices, obj.normals, obj.faces, ...
                                  obj.composition);
        end % End copy method
        
        function moveBy(obj, x)
        % Moves the TriagSurface object by the specidfed amount.
        %
        % INPUT:
        %  x - a 3 element row vector
            if ~isequal(size(x), [1 3])
                error('Can only move a sample by an amount [x y z]')
            end
            obj.vertices = obj.vertices + repmat(x, obj.nVertices, 1);
        end % End move function
        
        function plate_align(obj)
        % Correcly aligns the pinhole plate for the simulation. Should only be
        % used for that purpose.
            V = obj.vertices;
            N = obj.normals;
            
            obj.vertices(:,2) = -V(:,3);
            obj.vertices(:,3) = -V(:,2);
            obj.normals(:,2) = -N(:,3);
            obj.normals(:,3) = -N(:,2);
            
            obj.moveBy([0, max(abs(obj.vertices(:,2))), 0]);
        end % End pinhole plate align function
        
        function rotateY(obj)
        % Rotates the object by 90deg clockwise about the y axis.
            V = obj.vertices;
            N = obj.normals;
            
            obj.vertices(:,3) = V(:,1);
            obj.vertices(:,1) = -V(:,3);
            obj.normals(:,3) = N(:,1);
            obj.normals(:,1) = -N(:,3);
        end % End rotation function
        
        function rotateX(obj)
        % Rotates the object by 90deg clockwise about the x axis.
            V = obj.vertices;
            N = obj.normals;
            
            obj.vertices(:,2) = V(:,3);
            obj.vertices(:,3) = -V(:,2);
            obj.normals(:,2) = N(:,3);
            obj.normals(:,3) = -N(:,2);    
        end % End rotation function
        
        function rotateZ(obj)
        % Rotates the object by 90deg clockwise about the x axis.
            V = obj.vertices;
            N = obj.normals;
            
            obj.vertices(:,1) = V(:,2);
            obj.vertices(:,2) = -V(:,1);
            obj.normals(:,1) = N(:,2);
            obj.normals(:,2) = -N(:,1);    
        end % End rotation function
        
        function rotate(obj, rot_axis, angle)
        % Rotates the object counterclockwise by the given angle.
        % Doesn't work!!!!! (issue with the normals)
        %
        % INPUTS:
        %  rot_axis - string 'x', 'y', 'z', the axis that is to be rotated
        %             around
        %  angle    - angle, in degrees, that the object is to be rotated by
            V = obj.vertices;
            N = obj.vertices;
            V2 = zeros(size(V));
            N2 = zeros(size(N));
            rotMat = [cosd(angle), -sind(angle);
                      sind(angle),  cosd(angle)];        
            
            switch rot_axis
                case 'x'
                    V2(:,1) = V(:,1);
                    V2(:,2:3) = (rotMat*V(:,2:3)')';
                    N2(:,1) = N(:,1);
                    N2(:,2:3) = (rotMat*N(:,2:3)')';
                case 'y'
                    V2(:,2) = V(:,2);
                    V2(:,[1,3]) = (rotMat*V(:,[1,3])')';
                    N2(:,2) = N(:,2);
                    N2(:,[1,3]) = (rotMat*N(:,[1,3])')';
                case 'z'
                    V2(:,3) = V(:,3);
                    V2(:,1:2) = (rotMat*V(:,1:2)')';
                    N2(:,3) = N(:,3);
                    N2(:,1:2) = (rotMat*N(:,1:2)')';
                otherwise
                    error('Please specify a correct axis for rotation, x,y,z');
            end
            
            obj.vertices = V2;
            obj.normals = N2;
        end
        
        function reflectNormals(obj)
        % Reflects the normals in the object:
        % obj.n <- -obj.n
            obj.normals = -obj.normals;
        end
        
        function patchPlot(obj, new_fig, fname)
        % Produces a plot of the sample surface, if a file name is provided
        % then the plot is saved to that file. Can either plot in a new figure
        % or plot over an already present patch.
        % 
        % INPUTS:
        %  new_fig - bool, true => make a new figure, false => overlay on last
        %            used. If not provided defaults to true.
        %  fname   - file to save the figure to, if this is not provided then
        %            the figure is not saved
            if nargin == 1
                new_fig = true;
            end
            
            if new_fig
                figure
            end
            
            F = obj.faces;
            V = obj.vertices;
            
            patch('faces', F, 'vertices', V, 'FaceColor', [0.8 0.8 1.0], ...
                  'EdgeColor',       'black',        ...
                  'FaceLighting',    'gouraud',     ...
                  'AmbientStrength', 0.15);
            
            % Add a camera light, and tone down the specular highlighting
            if ~new_fig
                camlight('headlight');
            end
            material('dull');
            
            % Fix the axes scaling, and set a nice view angle
            axis('image');
            xlabel('x/mm')
            ylabel('y/mm')
            zlabel('z/mm')
            view([-5 5 5]);
            
            if nargin > 2
                saveas(gcf, fname, 'epsc');
            end
        end % End plotting function.
        
        function delete(obj)
            delete(obj);
        end % End delete function
    end % End methods
    
end % End classdef

