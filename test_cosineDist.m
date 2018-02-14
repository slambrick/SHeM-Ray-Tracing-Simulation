% Copyright (c) 2018, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
%
%
% This is a simple script that tests the sampling of a cosine distribution that
% is used in the SHeM simulation. It can easliy be modified, along with the mex
% file that it calls.

%% Parameters

clear

% The number of points to plot
N = 100000;

% Should results be saved to a text file
save_text = false;

% Name to save the results to 
fname = 'testingCosineDist.csv';

% Do we need to recompile the mex files
recompile = false;

% The normal about which to sample the ditribution
normal = [0 0 1];

%% Generate points

% Compile the mex file, comment out if this is not required
if recompile
    mex -lgsl -lgslcblas mexFiles/cosineMex.c mexFiles/small_functions.c %#ok<*UNRCH>
end

% Get the positions
positions = cosineMex(N, normal);
positions = positions';

%% Save a plot results

% Save the results to a file
if save_text
    fid = fopen(fname, 'w');
    
    for i_=1:N
        fprintf(fid, '%f,', positions(i_,:));
        fprintf(fid, '\n');
    end
    
    fclose(fid);
end

% Produce a 3D scatter plot of the results
figure(1)
scatter3(positions(:,1), positions(:,2), positions(:,3), '.')
axis equal
xlabel('x')
ylabel('y')
zlabel('z')



