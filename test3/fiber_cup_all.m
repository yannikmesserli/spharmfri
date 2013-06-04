% In this file, we are extracting the samples from the fiber cup data
% http://www.lnao.fr/spip.php?rubrique79

clear all;


% We will 
addpath('../common');
addpath('../nifti');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIRECTIONS
% Read the file containing the direction
directions = dlmread('../fiber_cup/diffusion_directions_corrected.txt');
% Get the sampled position:
[theta_tot, phi_tot, r] = cart2sph(directions(:,1), directions(:,2), directions(:,3));

% Take only one part:
phi = phi_tot(2:65);
theta = theta_tot(2:65);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRAINING PART
% Load the nii file
nii = load_nii('../fiber_cup/dwi-b1500.nii');

% some voxel positions, see http://www.lnao.fr/spip.php?article109
voxel_pos = {[51, 23, 1], [46, 21, 1], [51, 34, 1], [47, 32, 1], [41, 33, 1], [46, 38, 1], [44, 46, 1], [38, 48, 1], [31, 39, 1], [21, 48, 1], [17, 45, 1], [12, 40, 1], [11, 25, 1], [12, 17, 1], [24, 24, 1], [36, 9, 1]};

% Get all the sample for the voxel 1
t = repmat(voxel_pos{1}, 64, 1);
t = [t (2:65)' ];
% Get index of the ground:
base_i = sub2ind(size(nii.img), voxel_pos{1}(1), voxel_pos{1}(2), voxel_pos{1}(3), 1);
% Get index of all the different samples
index = sub2ind(size(nii.img), t(:,1), t(:,2), t(:,3), t(:,4));
% Consctruct the sample:
sample_test = ones( size(index) ) ./ log(  double(nii.img(index)) - double(nii.img(base_i)) );
% normalize it:
sample_test = (sample_test - min(sample_test));
sample_test = sample_test ./ max(sample_test);

% Construct the true diffusion direction:
fri.Locations = [ 0 0 ];
fri.Weights = [ 1 ];

% Then solve it:
h = kernelTrain(sample_test, fri, phi', theta');

% Save it in a mat file:
cart_cord1 = zeros(numel(voxel_pos), 3);
cart_cord1(1, : ) = sph2cart( fri.Locations(:, 1)', fri.Locations(:, 2)', fri.Weights);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST PART:
% same acquisition trials

for i = 2:numel(voxel_pos)

	t = repmat(voxel_pos{i}, 64, 1);
	t = [t (2:65)' ];

	index = sub2ind(size(nii.img), t(:,1), t(:,2), t(:,3), t(:,4));
	s = ones( size(index) ) ./ log(  double(nii.img(index)) - double(nii.img(base_i)) );
	% normalize it:
	s = (s - min(s));
	s = s ./ max(s);

	% Then solve it:
	fri_est = solveFRI(s, 1, phi', theta', h);

	cart_cord1(i, : ) = sph2cart( fri_est.Locations(:, 1)', fri_est.Locations(:, 2)', fri_est.Weights);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECOND SET:


t2 = repmat(voxel_pos{1}, 64, 1);
t2 = [t2 (67:130)' ];
base_i2 = sub2ind(size(nii.img), voxel_pos{1}(1), voxel_pos{1}(2), voxel_pos{1}(3), 66);
index2 = sub2ind(size(nii.img), t2(:,1), t2(:,2), t2(:,3), t2(:,4));
sample_test2 = ones( size(index2) ) ./ log(  double(nii.img(index2)) - double(nii.img(base_i2)) );
% normalize it:
sample_test2 = (sample_test2 - min(sample_test2));
sample_test2 = sample_test2 ./ max(sample_test2);

% Then solve it:
h2 = kernelTrain(sample_test2, fri, phi', theta');


% Save it in a mat file:
cart_cord2 = zeros(numel(voxel_pos), 3);
cart_cord2(1, : ) = sph2cart( fri.Locations(:, 1)', fri.Locations(:, 2)', fri.Weights);

for i = 2:numel(voxel_pos)

	t = repmat(voxel_pos{i}, 64, 1);
	t = [t (2:65)' ];

	index = sub2ind(size(nii.img), t(:,1), t(:,2), t(:,3), t(:,4));
	s = ones( size(index) ) ./ log(  double(nii.img(index)) - double(nii.img(base_i)) );
	% normalize it:
	s = (s - min(s));
	s = s ./ max(s);

	% Then solve it:
	fri_est = solveFRI(s, 1, phi', theta', h);

	cart_cord2(i, : ) = sph2cart( fri_est.Locations(:, 1)', fri_est.Locations(:, 2)', fri_est.Weights);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK TRUE POSITION: COMPARISON

save('cart_position.mat', 'cart_cord1', 'cart_cord2');

