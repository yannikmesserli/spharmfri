close all;
clear all;

addpath('./spharm');

% Parameters:
K = 2;% the number of dirac
N = 2*K; % the number of moments
M = 2*K; % Number of samples
check = 0; % Print some test?


    
% Our signal f(theta, varphi) = sum a_k delta(theta - theta_k, varphi - varphi_k ):
% theta is in [ 0  pi ]
% varphi is in [0 2 * pi ]

fri.Locations = [ sort( rand(1, K) * pi/2 )' sort(rand(1, K) * 2 * pi   )'];
fri.Weights = sort(rand(1, K) );





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTING \hat f_l^m,  Eq. (11):
%


% Compute some common variable:

u = exp(- 1i  .*  fri.Locations(:, 2)  );
s = sin( fri.Locations(:, 1)  );


% The spherical harmonics of degree l and order m:
f = zeros(N, N);
% The legendre polynomials for each degree l and each order m
P = zeros(N, N, K);


for n = 0:N-1
    Ptmp = legendre(n, cos(fri.Locations(:,1)) ); % normilzation to avoid the N coeff.
    P(n+1,1:n+1, 1:K) = Ptmp;
    if n == 0
        Ptmp = Ptmp';
    end
    
    for m = 0:n;
        f(n+1,m+1) = sum( fri.Weights .* (u.').^m .* Ptmp(m+1, :), 2);
    end 
    %f(n+1,1:n+1) = sum( repmat(fri.Weights, n+1, 1) .* repmat(u, 1, n+1).^(repmat( 0:n, K, 1)) .* Ptmp, 2);
end



% First we need only the diagonal: spherical of equal order and degree:
fnn = diag(f).';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT THE DIRAC:

%plot_sphere(K, fri);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPTUE THE SAMPLES:
s_test = samples(N, f, 0);

 

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
rts = roots([1 A.']);

% We can find phi:
phi_est     = sort( mod( - angle( rts), 2*pi));


print_check('Position check on phi: ',  phi_est, fri.Locations(:, 2));

% Theta:
theta_est_sign   =  sort(mod( real( asin(abs(rts))), pi/2));
print_check('Position test on theta: ',  theta_est_sign, fri.Locations(:, 1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND THE AMPLITUDE:
%

% Create the Vandermonde matrix:
V = repmat( rts, 1, K).'  .^ repmat( (0:K-1).', 1, K);

% Solve it:
a = sort( real(V \ z(1:K).'));


print_check('Amplitude check: ',  a, fri.Weights);







