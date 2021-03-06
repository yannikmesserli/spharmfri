% In this file I'll create a matrix reduced to only a full-rank matrix with the low-pass kernel. s
clear all;


% We will compute the samples with kernelP and f, fNeg
addpath('./spharm');


% Parameters:
K = 3;% the number of dirac
N = 2*K; % the number of moments

% Construct the signal:
fri.Locations = [  sort( rand(1, K) * pi )' sort(rand(1, K) * 2 * pi   )'];
fri.Weights = sort(rand(1, K) );

% Construct the matrix P
P = kernelP2(N);



% Compute the spherical harmonics:
[ftmp fNeg] = coeffFromFRI(fri);
% Get the diagonal only:
fnn = diag(ftmp);


% stest = samples(N, ftmp, 0);

s = P * fnn;

ftest = P \ s;

% Test:
if abs(ftest - fnn) > 0.00000001
	fprintf('Problem with the spherical harmonics \n');
end



