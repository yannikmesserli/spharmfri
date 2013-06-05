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

% Get all the sample for the voxel j
j = 3;
t = repmat(voxel_pos{j}, 64, 1);
t = [t (2:65)' ];
% Get index of the ground:
base_i = sub2ind(size(nii.img), voxel_pos{j}(1), voxel_pos{j}(2), voxel_pos{j}(3), 1);
% Get index of all the different samples
index = sub2ind(size(nii.img), t(:,1), t(:,2), t(:,3), t(:,4));
% Consctruct the sample:
s =  log( double(nii.img(base_i)) ) - log(  double(nii.img(index))  ) ;

% Construct the true diffusion direction:
load('../fiber_cup/tensor_feat.mat');
x = tensor_direction(j, :);
[p t w] = cart2sph(x(1), x(2), x(3));
fri.Locations = [ p t ];
fri.Weights = [ w ];

% Then solve it:
h = kernelTrain(sample_test, fri, phi', theta');

% Save it in a mat file:
cart_cord1 = zeros(length(voxel_pos), 3);
cart_cord1(j, : ) = x;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST PART:
% same acquisition trials

for i = 1:length(voxel_pos)
	if i ~= j
		fprintf('.');
		t = repmat(voxel_pos{i}, 64, 1);
		t = [t (2:65)' ];

		index = sub2ind(size(nii.img), t(:,1), t(:,2), t(:,3), t(:,4));
		s =  log( double(nii.img(base_i)) ) - log(  double(nii.img(index))  ) ;

		% Then solve it:
		fri_est = solveFRI(s, 1, phi', theta', h);

		[x y z] = sph2cart( fri_est.Locations(:, 1)', fri_est.Locations(:, 2)', fri_est.Weights);
		cart_cord1(i, : ) = [x y z];
	end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECOND SET:


t2 = repmat(voxel_pos{j}, 64, 1);
t2 = [t2 (67:130)' ];
base_i2 = sub2ind(size(nii.img), voxel_pos{j}(1), voxel_pos{j}(2), voxel_pos{j}(3), 66);
index2 = sub2ind(size(nii.img), t2(:,1), t2(:,2), t2(:,3), t2(:,4));
s2 =  log( double(nii.img(base_i2)) ) - log(  double(nii.img(index2))  ) ;

% Then solve it:
h2 = kernelTrain(s2, fri, phi', theta');


% Save it in a mat file:
cart_cord2 = zeros(length(voxel_pos), 3);
cart_cord2(j, : ) = x;

fprintf('\n');
for i = 1:length(voxel_pos)
	if i ~= j
		fprintf('.');
		t = repmat(voxel_pos{i}, 64, 1);
		t = [t (2:65)' ];

		index = sub2ind(size(nii.img), t(:,1), t(:,2), t(:,3), t(:,4));
		s = ones( size(index) ) ./ log(  double(nii.img(index)) - double(nii.img(base_i)) );
		% normalize it:
		s = (s - min(s));
		s = s ./ max(s);

		% Then solve it:
		fri_est = solveFRI(s, 1, phi', theta', h);

		[x y z] = sph2cart( fri_est.Locations(:, 1)', fri_est.Locations(:, 2)', fri_est.Weights);
		cart_cord2(i, : ) = [x y z];
	end

end
fprintf('\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK TRUE POSITION: COMPARISON
datetime=datestr(now);
datetime=strrep(datetime,':','_');
datetime=strrep(datetime,'-','_');
datetime=strrep(datetime,' ','_');

save( ['data/fiber_cup_all_' datetime], 'cart_cord1', 'cart_cord2');

