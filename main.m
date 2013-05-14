close all;
clear all;

addpath('./spharm');




% Create a random fri signal with 5 picks in [0 1] and we coefficient
% ranomly choosen between [-1 1]

% Parameters:
K = 4;% the number of dirac
N = 2*K; % the number of moments
M = 2*K; % Number of samples
check = 0; % Print some test?

    
% our signal f(theta, varphi) = sum a_k delta(theta - theta_k, varphi - varphi_k ):
% theta is in [ 0  pi ]
% varphi is in [0 2 * pi ]

fri.Locations = [ sort( rand(1, K) * pi )' sort(rand(1, K) * 2 * pi   )'];
fri.Weights = sort(rand(1, K) );


% Compute some common variable:

u = exp(- 1i  .*  fri.Locations(:, 2)  );
s = sin( fri.Locations(:, 1)  );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTING \hat f_l^m,  Eq. (11):
%

% The spherical harmonics of degree l and order m:
f = zeros(N, N);
% The legendre polynomials for each degree l and each order m
P = zeros(N, N, K);
% Normalizing factor:
N_fact = ones(N, N);


for n = 0:N-1
    Ptmp = legendre(n, cos(fri.Locations(:,1)) ); % normilzation to avoid the N coeff.
    P(n+1,1:n+1, 1:K) = Ptmp;
    if n == 0
        Ptmp = Ptmp';
    end
    % N_fact = 1;
    f(n+1,1:n+1) = sum( repmat(fri.Weights, n+1, 1) .* repmat( u.^n, 1, n+1).' .* Ptmp, 2);
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND PHI:
% Annihilating filter:
%


% First we need only the diagonal: spherical of equal order and degree:
fnn = diag(f).';
N_fact_nn = diag(N_fact).';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT THE DIRAC:

%plot_sphere(K, fri);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPTUE THE SAMPLES:
s_test = samples(N, f, 0);

 




% WHAT IS THIS NORMALIZATION FACTOR N??? -> use trick...

% Computing the true z from Eq. (13):
z_true = zeros(1, N);
for n = 0:N-1
    z_true(n+1) = sum(  fri.Weights' .* (s.* u).^n);
end

% Compute on each spherical harmonics what is the coefficient:
N_fact_real = real(fnn)./(real(z_true));
for n = 0:N-1
    % Normalizing factor of the Legendre polynomial N_n with repest to Rodriguez's fromula:
    t = @(n) (-1)^n * factorial(2*n-2)/factorial(n-1)*(n - 0.5)/2^(n-2);

    N_fact_real(n+1);
    if n > 0
        fprintf('%f = %f \n', t(n), N_fact_real(n+1) )
    else
        fprintf('%f = %f \n', 1.0, N_fact_real(n+1) )
    end
    
    %( (-1)^n * (n + 0.5) / factorial(2*n) )

    % Normalizing factor in physical convention (c.f. wikipedia)
    %sqrt( (2.*n + 1 )/ (4*pi * factorial(2.*n)  )  )
    % With the legendre normalizing coef:
    % sqrt( (2.*n + 1 )/ (4*pi * factorial(2.*n)  )  ) * (-1)^n * sqrt( (n + 0.5) / factorial(2*n) )
end
z = fnn ./ N_fact_real;
z = z + randn(size(z)).*0.2 ;


X = toeplitz( z(K:N-1), fliplr(  z(1:K) )  );
Y = - z(K+1:N).';
A = X \ Y;
rts = roots([1 A.']);

% We can find phi:
phi_est     = sort( mod( - angle( rts), 2*pi));


print_check('Position check on phi: ',  phi_est, fri.Locations(:, 2));


theta_est_sign   =  sort(mod( real( asin(abs(rts))), pi));

print_check('Position test on theta: ',  theta_est_sign, fri.Locations(:, 1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND THE AMPLITUDE:
%

% Create the Vandermonde matrix:
V = repmat( rts, 1, K).'  .^ repmat( (0:K-1).', 1, K);

% Solve it:
a = sort( real(V \ z(1:K).'));


print_check('Amplitude check: ',  a, fri.Weights);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND THETA:
% Annihilating filter:
%

% Compute w from Eq. (14):
w = zeros(K, 1);
c = cos( fri.Locations(:, 1) );
for n = 1:K
    w(n) = sum(  fri.Weights' .*  c .* (s.* u).^(n-1) );
end

% Create the Vandermonde matrix:
W = repmat( rts, 1, K).'  .^ repmat( (0:K-1)', 1, K) .* repmat( a', K, 1);

% Solve it:
c_est = (W \ w);


theta_est   =  sort(mod( real(acos(c_est)), pi));

print_check('Position check on theta: ', theta_est, fri.Locations(:, 1));




