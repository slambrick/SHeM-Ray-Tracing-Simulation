% Two sets of data
load('/home/sam/SMF_local/0073_deep_trench_test/scatteringData.mat')
deep_trench_data = simulationData;
load('/home/sam/SMF_local/0074_shallow_trench_test/scatteringData.mat')
shallow_trench_data = simulationData;

% Produce single scattering images with the same greyscale
s1 = reshape(deep_trench_data.counters{1}(1,:,:), [131, 50]);
s2 = reshape(shallow_trench_data.counters{1}(1,:,:), [131, 50]);
a1 = deep_trench_data.cntrSum{1};
a2 = shallow_trench_data.cntrSum{1};
norm_value = mean([a1(:,1); a2(:,1)]);
norm_fact1 = mean(s1(:,1))/norm_value;
norm_fact2 = mean(s2(:,1))/norm_value;
s1 = s1/norm_fact1;
s2 = s2/norm_fact2;

% Max/min intensity values
m1 = max([deep_trench_data.cntrSum{1}(:); shallow_trench_data.cntrSum{1}(:)]);
m2 = min([shallow_trench_data.cntrSum{1}(:); shallow_trench_data.cntrSum{1}(:)]);
% Lowest non-zero value in the single scattering
s_min = min([min(s1(s1 > 0)); min(s2(s2 > 0))]);
m_scale = [m2/2, m1];

im1 = mat2gray(s1, m_scale);
im2 = mat2gray(s2, m_scale);

% Plot images
im_a1 = deep_trench_data.imageAll('scale', 'manual', 'specifyScale', m_scale);
im_a2 = shallow_trench_data.imageAll('scale', 'manual', 'specifyScale', m_scale);
figure;
imshow(im1);
figure;
imshow(im2);

% Write images out to files
imwrite(im_a1, 'deep_multipleScattering.png')
imwrite(im_a2, 'shallow_multipleScattering.png')
imwrite(im1, 'deep_singleScattering.png')
imwrite(im2, 'shallow_singleScattering.png')
