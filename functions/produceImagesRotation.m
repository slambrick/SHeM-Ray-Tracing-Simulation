% produceImagesRotation.m
%
% Produces images of the sample for all the rotation cases.
function produceImagesRotation(simData, rot_angles, dataPath)
    for i_=1:length(simData)
        save_path = [dataPath '/rot' num2str(rot_angles(i_))];
        if ~exist(save_path, 'dir')
            mkdir(save_path);
        end
        simData{i_}.produceImages(save_path);
    end
end

