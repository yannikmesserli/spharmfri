

clear all;
close all;

% We will compute the samples with kernelP and f, fNeg
addpath('../common');

% Parameters:
K = 2;% the number of dirac
N = 2*K; % the number of moments

% Construct the signal:
fri.Locations = [  sort( rand(1, K) * pi )' sort(rand(1, K) * 2 * pi   )'];
fri.Weights = sort(rand(1, K) );



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

abs_diff = zeros(10, 1);
abs_diff_2 = zeros(10, 1);
hold on;
for power = 7:10
	for n = 1:10

		UP = n*16;
		phi 	= linspace(0, 2*pi - 2*pi/UP  , UP);  % Azi
		theta  	= linspace(0,   pi - pi/UP + 0.001   , UP);
		P = kernelP(N, phi, theta);

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


		s_noisy = s + 10^(-power) .* randn(size(s));


		f_noisy =  pinv(P) * s_noisy;

		abs_diff(n) = mean( abs(  (f_noisy - f_true)   ) );
		abs_diff_s(n) = mean( abs(  (s_noisy - s)   ) );

		

	end

	plot(1:10, abs_diff);
	
end
% legend(sprintf('Noise poser 10^%d', 5), sprintf('Noise poser 10^%d', 5), sprintf('Noise poser 10^%d', 5), sprintf('Noise poser 10^%d', 5), sprintf('Noise poser 10^%d', 5));
hold off;