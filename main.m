close all;
clear all;

addpath('./spharm');
addpath('./common');

% Parameters:
K = 2;% the number of dirac
N = 2*K; % the number of moments
check = 0; % Print some test?


    
% Our signal f(theta, varphi) = sum a_k delta(theta - theta_k, varphi - varphi_k ):
% theta is in [ 0  pi ]
% varphi is in [0 2 * pi ]

fri.Locations = [  [ 0.2 * pi  0.8 * pi ]' sort(rand(1, K) * 2 * pi   )'];
fri.Weights = sort(rand(1, K) );



% Compute the spherical coefficient for later use for this signal:
f = coeffFromFRI(fri);

% We'll need the diagonal: spherical of equal order and degree:
fnn = diag(f).';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT THE DIRAC:

%plot_sphere(K, fri);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPTUE THE SAMPLES:
s_test = samples(N, f, 0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE THE SPHERICAL COEFFICIENTS FROM THE SAMPLES:
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND PHI:
%

Nn = Normalizing(N);
% Construct z:
z = fnn ./ Nn;


% Annihilating filter:
X = toeplitz( z(K:N-1), fliplr(  z(1:K) )  );
Y = - z(K+1:N).';
A = X \ Y;
r = roots([1 A.']);

% We can find phi:
[phi_est, IX] = sort( mod( - angle(r), 2*pi));
print_check('Position check on phi: ',  phi_est, fri.Locations(:, 2));


% Put the r_k vector in the right order:
r = r(IX);
% Check if we have the correct order:
r_test(fri, (r));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND THE AMPLITUDE:
%

% Create the Vandermonde matrix:
V = repmat( r, 1, K).'  .^ repmat( (0:K-1).', 1, K);

% Solve it:
a = sort( real(V \ z(1:K).'));


print_check('Amplitude check: ',  a', fri.Weights);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND THETA:
% Annihilating filter:
%

% Compute w from Eq. (14):
w = zeros(K, 1);

for n = 1:K
    w(n) = -1 .* f(n+1, n) ./  Nn(n+1);
end

% Check if the vector w is correct:
w_test(fri, w);


% Create the Vandermonde matrix:
W = repmat( r.', K, 1)  .^ repmat( (0:K-1)', 1, K) .* repmat( a', K, 1);


% Solve it:
c_est = (W \ w);

% We can find theta:
theta_est   =  sort(mod( real(acos(c_est)), pi));

print_check('Position check on theta: ', theta_est, fri.Locations(:, 1));

