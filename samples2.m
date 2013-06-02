% We are using a low pass-filter 

clear all;


% We will compute the samples with kernelP and f, fNeg
addpath('./spharm');
% Parameters:
K = 1;% the number of dirac
N = 2*K; % the number of moments

% Construct the signal:
fri.Locations = [  [ 0.2 * pi  0.8 * pi ] sort(rand(1, K) * 2 * pi   )];
fri.Weights = sort(rand(1, K) );

% Construct the matrix P
P = kernelP(N);


f = zeros(length(P), 1);
% Compute the spherical harmonics:
[ftmp fNeg] = coeffFromFRI(fri);
% then reconstruct the vector in the right order:

% Dummy var to know where we are:
current = 0;
for l = 0:N-1
	
	current = current + max( 2*(l - 1)+1, 0);
	nbrep = 2*l + 1;
	
	f(current+1:current + (nbrep-1)/2 ) = fliplr(fNeg(l+1, 2:(nbrep-1)/2+1 ));
	f(current + (nbrep+1)/2 :current + nbrep ) = ftmp(l+1, 1:(nbrep+1)/2 );
end

% stest = samples(N, ftmp, 0);
rank(P)
cond(P)

s = P * f;

ftest = pinv(P) * s;


if abs(ftest - f) < 0.00000001
	fprintf('ok \n');
else
	fprintf('pas ok \n');
	ftest
	f
	sum(abs(ftest - f))
end