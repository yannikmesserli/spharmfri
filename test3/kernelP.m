% We are computing the the matrix P using Eq. (9) in the paper
% For arbitrary sampling kernel value h:

function P = kernelP(nb_dirac, PHI, THETA, h) 
	addpath('../spharm');

	L = 2 * nb_dirac;
	% The number of spherical coefficient:
	M = L^2;
	% Number of samples:
	N = length(PHI);
	
	if N >= M
		% Define our matrix P of size MxM
		P = zeros(M,N);

		% Dummy var to know where we are:
		current = 0;
		% fill the alternative matrix with P_l
		for l = 0:L-1
			Lmn 	= legendre(l, cos(THETA));
			for i = -l:1:l
				current = current +1;
				P(current, :) = constructP(l, i, THETA, PHI, h, Lmn, L);
				% fprintf('constructiong l=%d i=%d, \tY = %g %gi,  \tL = %f, \tNeg = %f\n', l, i, P(current, 1),  P(current, 1),  Lmn( abs(i)+1, 1), NSpher( abs(i), l));
			end
		end
		% Oups, make it right:
		P = P.';
	else
		throw( MException('kernelP:invalide', 'Not enough sample position'));
	end
end

% Compute the value of \hat p_l^m from Eq. (8) of the paper:
function P = constructP(l, m, THETA, PHI, h, Lmn, L)


	% For negative polynomial, I'm using this relationship: 
	% http://en.wikipedia.org/wiki/Associated_Legendre_polynomials#Negative_m_and.2For_negative_.E2.84.93
	% So:
	Nneg = @(m, l) (-1)^m * factorial(l-m)/factorial(l+m);

	% For the coefficient, we are using:
	NSpher = @(m, l) sqrt( ((2*l+1)/(4*pi)) *  factorial(l-m)/factorial(l+m)  );

	% Start value:
	P = 0;
	% Compute sum:
	for i = -l:1:l
		% Refactor the indices so that they fit the matrix indices 0 => L^2+1, -L^2 => 1, L^2 => 2L^2+1
		im = i*m + L^2 + 1;
		if i < 0
			P = P + ...
				NSpher( abs(i), l) * Nneg( abs(i), l) .* Lmn( abs(i)+1, :) .* exp( 1i * i * PHI) ... % Spherical harmonic Y_l^m
				.* h(im, l+1); % Sampling kernel coefficient
		elseif i >= 0
			P = P + ...
				NSpher( abs(i), l)					  .* Lmn( abs(i)+1, :) .* exp( 1i * i * PHI) ... % Spherical harmonic Y_l^m
				.* h(im, l+1); % Sampling kernel coefficient
		end
	end
	

end