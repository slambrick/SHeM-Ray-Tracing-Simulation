% re_rotate_images.m
%
% Copyright (c) 2020, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
%
% Generates the raster pattern for the 2D scan
%
% Calling syntax: 
%
% INPUTS:
%
% OUTPUTS:
function ims = re_rotate_images(simulationData, rot_angles, thePath)
    ims = cell(size(rot_angles));
    
    for i_=1:length(simulationData)
        simData = simulationData{i_};
        theta = rot_angles(i_);
        
        subPath = [thePath '/rotation' num2str(rot_angles(i_))];
        
        % Plot only the complete images on a series of figures
        figure;
        
        for j_=1:simData.n_detector
            single = imrotate(simData.getSingle(j_), theta - 90);
            imwrite(mat2gray(single), [subPath '/single_rotated' num2str(j_) '.png']);
    
            multiple = imrotate(simData.getMultiple(j_), theta - 90);
            imwrite(mat2gray(multiple), [subPath '/multiple_rotated' num2str(j_) '.png']);

            effuse = imrotate(simData.counter_effusive{j_}, theta - 90);
            imwrite(mat2gray(effuse), [subPath '/effuse_rotated' num2str(j_) '.png']);

            all = imrotate(simData.cntrSum{j_}, theta - 90);
            imwrite(mat2gray(all), [subPath '/all_rotated' num2str(j_) '.png']);
        
            subplot(1, simData.n_detector, j_)
            imshow(mat2gray(all))
            title(['Detector ' num2str(j_)]);
            
            im.Single{j_} = single;
            im.Multiple{j_} = multiple;
            im.Effuse{j_} = effuse;
            im.All{j_} = all;
        end
        
        ims{i_} = im;
    end
end

