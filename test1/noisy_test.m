

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

rmse_f = zeros(length(db_axis), 1);
rmse_angle = zeros(length(db_axis), 1);
rmse_weight = zeros(length(db_axis), 1);
rmse_all = zeros(length(db_axis), 1);




for i = 1:length(db_axis);

	db = db_axis(i);
	counter = 0;
	fprintf('%d / %d: ', i, length(db_axis) ) ;

	for n = 1:500
		
		if mod(n, 50) == 0
			fprintf('.');
		end

		

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
			

			% Then solve it:
			fri_est = solveFRI( [snn_noisy' sn1n_noisy'] );
			
			% comparte the two siganls:
			[a l all] = RMSE_FRI(fri_est, fri);
			% Save the mean square difference:
			rmse_weight(i) = rmse_weight(i) + l;
			rmse_angle(i) = rmse_angle(i) + a;
			rmse_all(i) = rmse_all(i) + all;
			rmse_f(i) = rmse_f(i) + RMSE( [ fnn_noisy' fn1n_noisy'], [fnn_true' fn1n_true' ] );


		end

	end
	rmse_f(i) = rmse_f(i) / counter;
	rmse_angle(i) = rmse_angle(i) / counter;
	rmse_weight(i) = rmse_weight(i) / counter;

	fprintf('\n');
end



% Save data:
%Build a name for the matfile
datetime=datestr(now);
datetime=strrep(datetime,':','_');
datetime=strrep(datetime,'-','_');
datetime=strrep(datetime,' ','_');

save( ['rmse_' datetime], 'rmse_angle', 'rmse_weight', 'rmse_f', 'legends', 'db_axis');

