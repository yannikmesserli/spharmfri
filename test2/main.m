

clear all;


% We will compute the samples with kernelP and f, fNeg
addpath('../common');

% Parameters:
K = 2;% the number of dirac
N = 2*K; % the number of moments

% Construct the signal:
fri.Locations = [  sort( rand(1, K) * pi )' sort(rand(1, K) * 2 * pi   )'];
fri.Weights = sort(rand(1, K) );



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE THE SAMPLES:
%


% Construct the matrix P
% Which consist of only the independant column of the 
% kernel of a low-pass filter:
% UP = 40;
% phi 	= linspace(0, 2*pi - randn(1) * 2*pi/UP  , UP);  % Azi
% theta  	= linspace(0,   pi - randn(1) * pi/UP, UP);
% P = kernelP(K, phi, theta);
P = kernelP(K);


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

% power = 2;
% s_noisy = s + 10^(-power) .* randn(size(s)) + 10^(-power) .* randn(size(s)) * 1i;
% s_noisy

% f_noisy =  pinv(P) * s_noisy;

% abs_diff_f = mean( abs(  (f_noisy - f_true)   ) )
% abs_diff_s = mean( abs(  (s_noisy - s)   ) )

% Then solve it:
fri_est = solveFRI(s, 2);

print_check('Position check on phi: ',  fri_est.Locations(:, 2), fri.Locations(:, 2));
print_check('Position check on theta: ', fri_est.Locations(:, 1), fri.Locations(:, 1));
print_check('Amplitude check: ',  fri_est.Weights, fri.Weights);
