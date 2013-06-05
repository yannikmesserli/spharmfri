function [rmse_a rmse_l rmse_all] = RMSE_FRI(fri1, fri2)
	rmse_l = sqrt(mean( abs(fri1.Weights - fri2.Weights).^2 ));

	% Transform to cartesian value with length 1.0 to compare everything:
	[x y z] = sph2cart( fri1.Locations(:, 1)', fri1.Locations(:, 2)', ones(size(fri1.Locations(:, 2)')) );
	[x2 y2 z2] = sph2cart( fri2.Locations(:, 1)', fri2.Locations(:, 2)', ones(size(fri2.Locations(:, 2)')) );
	% Concatenate everything in a big vector:
	vect1 = [x y z];
	vect2 = [x2 y2 z2];
	rmse_a = sqrt(mean( abs(vect1 - vect2).^2 ));


	% Transform to cartesian value with length 1.0 to compare everything:
	[x y z] = sph2cart( fri1.Locations(:, 1)', fri1.Locations(:, 2)',  fri1.Weights);
	[x2 y2 z2] = sph2cart( fri2.Locations(:, 1)', fri2.Locations(:, 2)', fri2.Weights );
	% Concatenate everything in a big vector:
	vect1 = [x y z];
	vect2 = [x2 y2 z2];
	rmse_all = sqrt(mean( abs(vect1 - vect2).^2 ));
end