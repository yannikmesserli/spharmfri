

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

for j = 1:2
	UP = 16*j;
	counter = 0;

	for i = 1:length(db_axis);
		for n = 1:500
		
			db = db_axis(i);

			

			
			phi 	= linspace(0, 2*pi - 2*pi/UP  , UP); 
			theta  	= linspace(0,   pi - randn(1) * pi/UP, UP);
			P = kernelP(K, phi, theta);

			% Compute the spherical harmonics:
			[ftmp ftmpNeg] = coeffFromFRI(fri);
			% Pick up only the ones in the diagonal
			f_true = spharm2vect(ftmp, ftmpNeg);


			% Compute the samples:
			s = P * f_true;


			% Test to be sure our matrix is not blowing up the system:
			ftest = pinv(P) * s;
			if abs(ftest - f_true) > 0.00000001
				fprintf('Problem with the matrix inversion \n');
			end

			
			s_noisy = addNoise(s, db);
			

			f_noisy =  pinv(P) * s_noisy;
			counter = counter +1;
			rmse_f(j, i) = rmse_f(j, i) + RMSE(f_noisy, f_true);


			fri_est = solveFRI(s_noisy, K, phi, theta);
			
			% Reshape as a unique vector:
			t1 = reshape(fri_est.Locations, 1, numel(fri_est.Locations));
			t2 = reshape(fri.Locations, 1, numel(fri.Locations));

			% Save the mean square difference:
			rmse_l(j, i) = rmse_l(j, i) + RMSE( t1, t2);
			rmse_a(j, i) = rmse_a(j, i) + RMSE(fri_est.Weights, fri.Weights);

		end
		rmse_f(j, i) = rmse_f(j, i) / counter;
		rmse_a(j, i) = rmse_a(j, i) / counter;
		rmse_l(j, i) = rmse_l(j, i) / counter;
	end
	legends{j} = sprintf('%d x nb. min of samples', j);
end


% Save data:
%Build a name for the matfile
datetime=datestr(now);
datetime=strrep(datetime,':','_');
datetime=strrep(datetime,'-','_');
datetime=strrep(datetime,' ','_');

save( ['rmse_' datetime], 'rmse_a', 'rmse_l', 'rmse_f', 'legends', 'db_axis');

