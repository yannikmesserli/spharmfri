% COMPUTING \hat f_l^m,  Eq. (11):
function [f fneg P] = coeffFromFRI(fri)
    
    % Parameters:
    K = length(fri.Weights);% the number of dirac
    N = 2*K; % the number of moments

    % Compute some common variable:
    u = exp(- 1i  .*  fri.Locations(:, 2)  );
    s = sin( fri.Locations(:, 1)  );


    % The spherical harmonics of degree l and order positive m:
    f = zeros(N, N);
    % The spherical harmonics of degree l and order negative m:
    fneg = zeros(N, N);
    % The legendre polynomials for each degree l and each order m
    P = zeros(N, N, K);

    
    % For negative polynomial, I'm using this relationship: 
    % http://en.wikipedia.org/wiki/Associated_Legendre_polynomials#Negative_m_and.2For_negative_.E2.84.93
    % So:
    M = @(m, l) (-1)^m * factorial(l-m)/factorial(l+m);

    for n = 0:N-1
        Ptmp = legendre(n, cos(fri.Locations(:,1)) ); % normilzation to avoid the N coeff.
        P(n+1,1:n+1, 1:K) = Ptmp;
        if n == 0
            Ptmp = Ptmp';
        end
        
        for m = 0:n;
            f(n+1,m+1) = sum( fri.Weights .* (u.').^m .* Ptmp(m+1, :), 2);
            fneg(n+1,m+1) = M(m, n) * sum( fri.Weights .* (u.').^m .* Ptmp(m+1, :), 2);
        end 
        %f(n+1,1:n+1) = sum( repmat(fri.Weights, n+1, 1) .* repmat(u, 1, n+1).^(repmat( 0:n, K, 1)) .* Ptmp, 2);
    end

end