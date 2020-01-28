% TriagSurface.m
%
% Copyright (c) 2018-19, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the
% GNU/GPL-3.0-or-later.
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
%  compositions- A length-m cell array containing the material name of the respective face
%  materials   - A map from material names to their properties:
%                color, scattering function and its parameters
%  nTriag      - The number of triangles in the surface mesh.
%  nVertices   - The number of veritces in the surface mesh.
classdef TriagSurface < handle

    properties (SetAccess = private)
        vertices    % Vertices in the triangulation
        faces       % Contains the reference the vertices that make up the
                    % triangles
        normals     % Normals to the triangles
        compositions% The material name for each triangle
        materials   % the material library holding descriptions of the materials
        nTriag      % The number of triangles in the surface
        nVertices   % The number of vertices in the surface
    end % End properties

    methods
        function obj = TriagSurface(vertices, fdef, fnorm, fmat, materials)
        % Constructor for TriagSurface class. Can be called with 0 or 4
        % arguments. 0 arguments creats an empty object.
        %
        % INPUTS:
        %  vertices - The coordinates of the vertices.
        %  fnorm - The normals to the faces of the triangles.
        %  fdef - The 'faces', contains the index of the vertices that make up
        %         each face.
        %  fmat - the material name of each face
        %  materials - the materials library
        %
        % OUTPUT:
        %  obj - A TriagSurface object that contains the surface.
            if (nargin == 0)
                % Default, creates an empty object
            elseif (nargin == 5)
                if size(size(vertices), 2) ~= 2
                    error('Input arguments must all be 2D arrays')
                end
                if (size(vertices,2) ~= 3 || size(fnorm,2) ~= 3 || size(fdef,2) ~= 3)
                    error('vertices, fnorm, fdef inputs must all have three columns')
                end
                if size(vertices,1) < 3
                    error('There must be at least three vertices in the surface')
                end
                if (size(fdef, 1) ~= size(fnorm,1) || size(fdef,1) ~= size(fmat, 1))
                    error('Inputs fdef, fnorm, fmat must have the same length')
                end
                if length(materials) < 1
                    error('There must be at least one defined material - the default')
                end
                obj.vertices = vertices;
                obj.faces = fdef;
                obj.normals = fnorm;
                obj.compositions = fmat;
                obj.materials = materials;
                obj.nTriag = size(fdef, 1);
                obj.nVertices = size(vertices, 1);
            else
                error('Wrong number of input arguments');
            end
        end % End constructor

        function newobj = copy(obj)
        % Copys the object and returns the copy.
            newobj = TriagSurface(obj.vertices, obj.faces, obj.normals, ...
                                  obj.compositions, obj.materials);
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

        function plate_align_newmicro(obj)
            % Offset correctly
            offset = min(obj.vertices(:,2));
            obj.moveBy([0, -offset, 0]);
        end

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

        function rotateGeneral(obj, axis, theta)
        % Arbitray rotation about any axis.
        %
        % axis  - string, 'x', 'y', 'z'
        % theta - double, rotation angle in degrees
           theta = theta*pi/180;
           s = sin(theta);
           c = cos(theta);
           switch axis
               case 'x'
                   R = [1, 0, 0; 0, c, -s, 0, s, c];
               case 'y'
                   R = [c, 0, s; 0, 1, 0; -s, 0, c];
               case 'z'
                   R = [c, -s, 0; s, c, 0; 0, 0, 1];
           end
           obj.vertices = (R*obj.vertices')';
           obj.normals = (R*obj.normals')';
        end

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

            unqiue_comps = unique(obj.compositions);
            for idx = 1:length(unqiue_comps)
                % group all faces with this composition, extract the color and plot.
                comp = unqiue_comps{idx};
                indexes = strcmp(obj.compositions, comp);
                F = obj.faces(indexes, :);
                C = obj.materials(comp); C = C.color;

                patch('faces', F, 'vertices', obj.vertices, 'FaceColor', C, ...
                  'EdgeColor',       'black',        ...
                  'FaceLighting',    'gouraud',     ...
                  'AmbientStrength', 0.15);
            end

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
