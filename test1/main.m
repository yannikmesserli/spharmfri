

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

% Test to be sure our matrix is not blowing up the system:
ftest = pinv(P) * snn;
if abs(ftest - fnn_true) > 0.00000001
	fprintf('Problem with the matrix inversion \n');
end


% Then solve it:
fri_est = solveFRI([snn' sn1n'] );

print_check('Position check on phi: ',  fri_est.Locations(:, 2), fri.Locations(:, 2));
print_check('Position check on theta: ', fri_est.Locations(:, 1), fri.Locations(:, 1));
print_check('Amplitude check: ',  fri_est.Weights, fri.Weights);
