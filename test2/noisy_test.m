

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
rmse_angle = zeros(2, length(db_axis));
rmse_weight = zeros(2, length(db_axis));
rmse_all = zeros(2, length(db_axis));
fri_all = zeros(2, length(db_axis), 500, 3);

for j = 1:2
	UP = 16*j;
	% Construct the sampling grid:
	phi 	= linspace(0, 2*pi - 2*pi/UP  , UP); 
	theta  	= linspace(0,   pi - randn(1) * pi/UP, UP);
	% Construct the matrix
	P = kernelP(K, phi, theta);

	% Compute the spherical harmonics:
	[ftmp ftmpNeg] = coeffFromFRI(fri);
	% Pick up only the ones in the diagonal
	f_true = spharm2vect(ftmp, ftmpNeg);

	for i = 1:length(db_axis);

		fprintf('%d / %d: ', i, length(db_axis) ) ;
		counter = 0;
		db = db_axis(i);

		for n = 1:500
		
			if mod(n, 50) == 0
				fprintf('.');
			end
			
			
			


			% Compute the samples:
			s = P * f_true;
			% Add noise
			s_noisy = addNoise(s, db);
			% Compute the sphericla back:
			f_noisy = pinv(P) * s_noisy;
			% Solve it
			fri_est = solveFRI(s_noisy, K, phi, theta);
			
			% comparte the two siganls:
			[a l all] = RMSE_FRI(fri_est, fri);
			% Save the mean square difference:
			rmse_weight(j, i) = rmse_weight(j, i) + l;
			rmse_angle(j, i) = rmse_angle(j, i) + a;
			rmse_all(j, i) = rmse_all(j, i) + all;
			rmse_f(j, i) = rmse_f(j, i) + RMSE(f_noisy, f_true);

		end

		rmse_f(j, i) = rmse_f(j, i) / counter;
		rmse_angle(j, i) = rmse_angle(j, i) / counter;
		rmse_weight(j, i) = rmse_weight(j, i) / counter;
		rmse_all(j, i) = rmse_all(j, i) / counter;

		fprintf('\n');
	end
	legends{j} = sprintf('%d x nb. min of samples', j);
end


% Save data:
%Build a name for the matfile
datetime=datestr(now);
datetime=strrep(datetime,':','_');
datetime=strrep(datetime,'-','_');
datetime=strrep(datetime,' ','_');

save( ['rmse_' datetime], 'rmse_angle', 'rmse_weight', 'rmse_f', 'legends', 'db_axis', 'fri_all');

