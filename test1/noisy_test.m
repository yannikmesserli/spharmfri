

clear all;
close all;

% We will compute the samples with kernelP and f, fNeg
addpath('../common');

% Parameters:
K = 2;% the number of dirac
N = 2*K; % the number of moments

% Construct the signal:
% fri.Locations = [  sort( rand(1, K) * pi )' sort(rand(1, K) * 2 * pi   )'];
fri.Locations = [ [0.2 * pi  0.8 * pi ]' [ 0.1 * 2 * pi 0.6 * 2 * pi ]'  ]; % Make well separated angle
fri.Weights = sort(rand(1, K) );



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

abs_diff = zeros(10, 1);
abs_diff_2 = zeros(10, 1);

legends = {};
db_axis = 10:0.2:80;

rmse_f = zeros(2, length(db_axis));
rmse_a = zeros(2, length(db_axis));
rmse_l = zeros(2, length(db_axis));





for i = 1:length(db_axis);

	counter = 0;
	fprintf('%d / %d: ', i, length(db_axis) ) ;

	for n = 1:500
		
		if mod(n, 50) == 0
			fprintf('.');
		end

		db = db_axis(i);

		% Construct the matrix P
		% Which consist of only the independant column of the 
		% kernel of a low-pass filter:
		P = kernelP(N);


		% Compute the spherical harmonics:
		[ftmp fNeg] = coeffFromFRI(fri);
		% Pick up only the ones in the diagonal
		fnn_true = diag(ftmp);
		% and the one just below:
		fn1n_true = diag(ftmp, -1);



		% Compute the samples:
		snn = P * fnn_true;
		% The first line is f_0^0 so remove it:
		sn1n = P(2:end, 2:end) * fn1n_true;

		
		% Check if it is ok
		fnn_est = pinv(P) * snn;
	

		if abs(fnn_est - fnn_true) < 0.00000001

			counter = counter + 1;

			% Add noise
			snn_noisy = addNoise(snn, db);
			sn1n_noisy = addNoise(sn1n, db);

			% Get back the coeff.
			fnn_noisy =  pinv(P) * snn_noisy;
			fn1n_noisy =  pinv(P(2:end, 2:end)) * sn1n_noisy;

			rmse_f(i) = rmse_f(i) + RMSE( [ fnn_noisy' fn1n_noisy'], [fnn_true' fn1n_true' ] );

			% Compute the error done:

			% Then solve it:
			fri_est = solveFRI( [snn_noisy' sn1n_noisy'] );
			
			% Reshape as a unique vector:
			t1 = reshape(fri_est.Locations, 1, numel(fri_est.Locations));
			t2 = reshape(fri.Locations, 1, numel(fri.Locations));

			% Save the mean square difference:
			rmse_l(i) = rmse_l(i) + RMSE( t1, t2);
			rmse_a(i) = rmse_a(i) + RMSE(fri_est.Weights, fri.Weights);

		end

	end
	rmse_f(i) = rmse_f(i) / counter;
	rmse_a(i) = rmse_a(i) / counter;
	rmse_l(i) = rmse_l(i) / counter;

	fprintf('\n');
end



% Save data:
%Build a name for the matfile
datetime=datestr(now);
datetime=strrep(datetime,':','_');
datetime=strrep(datetime,'-','_');
datetime=strrep(datetime,' ','_');

save( ['rmse_' datetime], 'rmse_a', 'rmse_l', 'rmse_f', 'legends', 'db_axis');

