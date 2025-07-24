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

    properties %(SetAccess = private)
        vertices    % Vertices in the triangulation
        faces       % Contains the reference the vertices that make up the
                    % triangles
        normals     % Normals to the triangles
        lattice     % Reciprocal lattice vectors
        compositions% The material name for each triangle
        materials   % the material library holding descriptions of the materials
        nTriag      % The number of triangles in the surface
        nVertices   % The number of vertices in the surface
    end % End properties

    methods
        function obj = TriagSurface(vertices, fdef, fnorm, lattice, fmat, materials)
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
            elseif (nargin == 6)
                if size(size(vertices), 2) ~= 2
                    error('Input arguments must all be 2D arrays')
                end
                if (size(vertices,2) ~= 3 || size(fnorm,2) ~= 3 || size(fdef,2) ~= 3)
                    error('vertices, fnorm, fdef inputs must all have three columns')
                end
                if (size(lattice,2) ~= 6)
                    error('lattice input must have 6 columns');
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
                obj.lattice = lattice;
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
            mat = containers.Map(...
                keys(obj.materials), ...
                values(obj.materials));
            newobj = TriagSurface(obj.vertices, obj.faces, obj.normals, obj.lattice, ...
                                  obj.compositions, mat);
        end % End copy method

        function moveBy(obj, x)
        % Moves the TriagSurface object by the specidfed amount.
        %
        % INPUT:
        %  x - a 3 element row vector
            if ~isequal(size(x), [1 3])
                error(['Can only move a sample by an amount [x y z]. ' x ' is invalid'])
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

        function rotateY(obj, ang)
        % Rotates the object by 90deg clockwise about the y axis.
            if nargin == 1
                ang = 90;
            end
            rotateGeneral(obj, 'y', ang);
        end % End rotation function

        function rotateX(obj, ang)
        % Rotates the object by 90deg clockwise about the x axis.
            if nargin == 1
                ang = 90;
            end
            rotateGeneral(obj, 'x', ang);
        end % End rotation function

        function rotateZ(obj, ang)
        % Rotates the object by 90deg clockwise about the z axis.
            if nargin == 1
                ang = 90;
            end
            rotateGeneral(obj, 'z', ang);
        end % End rotation function

        function rotateGeneral(obj, axis, theta, how)
        % Arbitray rotation about any axis.
        %
        % axis  - string, 'x', 'y', 'z'
        % theta - double, rotation angle in degrees
            if nargin == 3
                how = 'both';
            end
            theta = theta*pi/180;
            s = sin(theta);
            c = cos(theta);
            switch axis
                case 'x'
                    R = [1, 0, 0; 0, c, -s; 0, s, c];
                case 'y'
                    R = [c, 0, s; 0, 1, 0; -s, 0, c];
                case 'z'
                    R = [c, -s, 0; s, c, 0; 0, 0, 1];
            end
            if strcmp(how, 'both') || strcmp(how, 'geometry')
                obj.vertices = (R*obj.vertices')';
                obj.normals = (R*obj.normals')';
            end
            if strcmp(how, 'both') || strcmp(how, 'lattice')
                obj.lattice(:,1:3) = (R*obj.lattice(:,1:3)')';
                obj.lattice(:,4:6) = (R*obj.lattice(:,4:6)')';
            end
        end
        
        function rotateLattice(obj, axis, theta)
            rotateGeneral(obj, axis, theta, 'lattice');
        end

        function reflectNormals(obj)
        % Reflects the normals in the object:
        % obj.n <- -obj.n
            obj.normals = -obj.normals;
        end

        function reflect_axis(obj, axis_name)
            % reflect the coordinates on one axis
            % axis_name should be 'x', 'y' or 'z'
            switch axis_name
            case 'x'
                idx = 1;
            case 'y'
                idx = 2;
            case 'z'
                idx = 3;
            otherwise
                error('Axis name must be x, y or z');
            end

            V = obj.vertices;
            N = obj.normals;

            obj.vertices(:, idx) = -V(:, idx);
            obj.normals(:, idx) = -N(:, idx);
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
                C = obj.materials(comp);
                C = C.color;

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
        
        function normalise_lattice(obj) 
            for i_=1:obj.nTriag
                b1 = obj.lattice(i_,1:3);
                b2 = obj.lattice(i_,4:6);
                obj.lattice(i_,1:3) = b1/sqrt(sum(b1.^2));
                obj.lattice(i_,4:6) = b1/sqrt(sum(b2.^2));
            end
        end

        function delete(obj)
            delete(obj);
        end % End delete function
    end % End methods

end % End classdef
