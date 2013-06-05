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

figure; plotonsphere2(phi, theta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAMPLES
% Load the nii file
nii = load_nii('../fiber_cup/dwi-b1500.nii');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK TRUE POSITION
% Load the nii file
nii_tensor = load_nii('../fiber_cup/diff/dti_v1.nii');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% some voxel positions, see http://www.lnao.fr/spip.php?article109
voxel_pos = {[51, 23, 1], [46, 21, 1], [51, 34, 1], [47, 32, 1], [41, 33, 1], [46, 38, 1], [44, 46, 1], [38, 48, 1], [31, 39, 1], [21, 48, 1], [17, 45, 1], [12, 40, 1], [11, 25, 1], [12, 17, 1], [24, 24, 1], [36, 9, 1]};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



est_direction = zeros(length(voxel_pos), 3);
tensor_direction = zeros(length(voxel_pos), 3);
rmse = zeros(length(voxel_pos), 1);
fri_all = {};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST PART:
% same acquisition trials

for i = 1:numel(voxel_pos)
	fprintf('.');
	% Get all the sample for the voxel 1
	t = repmat(voxel_pos{i}, 64, 1);
	t = [t (2:65)' ];
	% Get index of the ground:
	base_i = sub2ind(size(nii.img), voxel_pos{i}(1), voxel_pos{i}(2), voxel_pos{i}(3), 1);
	% Get index of all the different samples
	index = sub2ind(size(nii.img), t(:,1), t(:,2), t(:,3), t(:,4));
	% Consctruct the sample:
	s =   log( double(nii.img(base_i)) ) - log(  double(nii.img(index))  ) ;

	% Then solve it:
	fri_est = solveFRI(s, 1, phi', theta');
	% Save it
	fri_all{i} = fri_est;
	% Save the estimate position in cartesian coord:
	est_direction(i, : ) = sph2cart( fri_est.Locations(:, 1)', fri_est.Locations(:, 2)', fri_est.Weights);
	% Save the estimate position with the tensor:
	tensor_direction(i, : ) = nii_tensor.img(voxel_pos{i}(1), voxel_pos{i}(2), voxel_pos{i}(3), 1);
	% Compute the error:
	rmse(i) = RMSE( est_direction(i, : ), tensor_direction(i, : ) );
end
fprintf('\n');

%Build a name for the matfile
datetime=datestr(now);
datetime=strrep(datetime,':','_');
datetime=strrep(datetime,'-','_');
datetime=strrep(datetime,' ','_');

save( ['fiber_cup_' datetime], 'rmse', 'est_direction', 'tensor_direction', 'fri_all');

