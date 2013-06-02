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