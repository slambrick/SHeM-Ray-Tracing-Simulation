function [I, I_avs, I_diff, I_comb_diff, f] = difference_image(simData)
    clc
    
    I = simData.cntrSum;
    max_max = max(max(cell2mat(I)));
    min_min = min(min(cell2mat(I)));
        
    % plot the 4 detectors with the same greyscale
    f{1} = figure;
    ind = [1,3,2,4];
    for i_=1:4
        subplot(2, 2, ind(i_))
        imshow(I{i_}, [min_min, max_max])
        title(['Detector ' num2str(i_)])
        disp(['Detector ' num2str(i_) ' std: ' num2str(std(I{i_}(:)))]);
    end
    print('Orignal_images.png', '-dpng', '-r600')
    
    % Plot the combined image
    f{2} = figure;
    I_av = zeros(size(I{1}));
    for i_=1:length(I)
        I_av = I_av + I{i_};
    end
    imshow(mat2gray(I_av));
    pause(0.05);
    title('Totally averaged image')
    disp(['Totla average std: ' num2str(std(I_av(:)))])
    print('Average_plot.png', '-dpng', '-r600')
    
    % Plot combined from 2 detectors next to each other
    f{3} = figure;
    I_avs{1} = I{1} + I{2};
    I_avs{2} = I{2} + I{4};
    I_avs{3} = I{1} + I{3};
    I_avs{4} = I{3} + I{4};
    max_max = max(max(cell2mat(I_avs)));
    min_min = min(min(cell2mat(I_avs)));
    
    subplot(2, 2, 1)
    imshow(I_avs{1}, [min_min, max_max])
    title('Detectors 1&2')
    disp(['Detectors 1&2 std: ' num2str(std(I_avs{1}(:)))]);
    subplot(2, 2, 2)
    imshow(I_avs{2}, [min_min, max_max])
    title('Detectors 2&4')
    disp(['Detectors 2&4 std: ' num2str(std(I_avs{2}(:)))]);
    subplot(2, 2, 3)
    imshow(I_avs{3}, [min_min, max_max])
    title('Detectors 1&3')
    disp(['Detectors 1&3 std: ' num2str(std(I_avs{3}(:)))]);
    subplot(2, 2, 4)
    imshow(I_avs{4}, [min_min, max_max])
    title('Detectors 3&4')
    disp(['Detectors 3&4 std: ' num2str(std(I_avs{4}(:)))]);
    print('Combined_plot.png', '-dpng', '-r600')
    
    % Plot combined from 2 detectors opposite from each other
    f{4} = figure;
    I_avs2{1} = I{1} + I{4};
    I_avs2{2} = I{2} + I{3};
    max_max = max(max(cell2mat(I_avs2)));
    min_min = min(min(cell2mat(I_avs2)));
    
    subplot(1, 2, 1)
    imshow(I_avs2{1}, [min_min, max_max])
    title('Detectors 1&4')
    disp(['Detectors 1&4 std: ' num2str(std(I_avs2{1}(:)))]);
    subplot(1, 2, 2)
    imshow(I_avs2{2}, [min_min, max_max])
    title('Detectors 2&3')
    disp(['Detectors 2&3 std: ' num2str(std(I_avs2{2}(:)))]);
    print('Combined_opposite_plot.png', '-dpng', '-r600')
    
    % Plot combining 3 detectors
    f{5} = figure;
    I_avs3{1} = I{1} + I{2} + I{3};
    I_avs3{2} = I{1} + I{2} + I{4};
    I_avs3{3} = I{2} + I{4} + I{3};
    I_avs3{4} = I{1} + I{3} + I{4};
    max_max = max(max(cell2mat(I_avs3)));
    min_min = min(min(cell2mat(I_avs3)));
    subplot(2, 2, 1)
    imshow(I_avs3{1}, [min_min, max_max])
    title('Detectors 2,1,3')
    disp(['Detectors 2,1,3 std: ' num2str(std(I_avs3{1}(:)))]);
    subplot(2, 2, 2)
    imshow(I_avs3{2}, [min_min, max_max])
    title('Detectors 4,2,1')
    disp(['Detectors 4,2,1 std: ' num2str(std(I_avs3{2}(:)))]);
    subplot(2, 2, 3)
    imshow(I_avs3{3}, [min_min, max_max])
    title('Detectors 3,4,2')
    disp(['Detectors 3,4,2 std: ' num2str(std(I_avs3{3}(:)))]);
    subplot(2, 2, 4)
    imshow(I_avs3{4}, [min_min, max_max])
    title('Detectors 1,3,4')
    disp(['Detectors 1,3,4 std: ' num2str(std(I_avs3{4}(:)))]);
    print('3Combined_plot.png', '-dpng', '-r600')
    
    % Plot two difference images from opposite angles
    f{6} = figure('Position', [0 0 700 400]);
    I_diff{1} = I{1} - I{2};
    subplot(1, 3, 1)
    imagesc(I_diff{1});
    colormap('parula')
    axis equal
    axis tight
    axis off
    title('Difference: 1 - 2')
    disp(['Difference 1 - 2 std: ' num2str(std(I_diff{1}(:)))]);
    I_diff{2} = I{1} - I{3};
    subplot(1, 3, 2)
    imagesc(I_diff{2});
    colormap('parula')
    axis equal
    axis tight
    axis off
    title('Difference: 1 - 3')
    disp(['Difference 1 - 3 std: ' num2str(std(I_diff{2}(:)))]);
    I_diff{3} = I{1} - I{4};
    subplot(1, 3, 3)
    imagesc(I_diff{3});
    colormap('parula')
    axis equal
    axis tight
    axis off
    title('Difference: 1 - 4')
    disp(['Difference 1 - 4 std: ' num2str(std(I_diff{3}(:)))]);
    print('Difference_plot.png', '-dpng', '-r600')
    
    colormap('gray')
    print('Difference_plot_grey.png', '-dpng', '-r600')
    
    % Plot the combined opposite difference images
    f{7} = figure;
    I_comb_diff{1} = I_avs{1} - I_avs{4};
    I_comb_diff{2} = I_avs{2} - I_avs{3};
    subplot(1, 2, 1)
    imshow(mat2gray(I_comb_diff{1}))
    disp(['Combined difference (1) std: ' num2str(std(I_comb_diff{1}(:)))]);
    subplot(1, 2, 2)
    imshow(mat2gray(I_comb_diff{2}))
    disp(['Combined difference (1) std: ' num2str(std(I_comb_diff{2}(:)))]);
    print('combined_difference.png', '-dpng', '-r600')
end

function save_SameScale(simData)
    I = simData.cntrSum;
    
    figure
    max_max = max(max(cell2mat(I)));
    min_min = min(min(cell2mat(I)));
    for i_=1:length(I)
        I2 = (I{i_} - min_min)/(max_max - min_min);
        imwrite(I2, [output_dir '/orignal_image' num2str(i_) '.png'], 'png')
    end
end