% We are computing the the matrix P using Eq. (8) in the paper. 
% - Our kernel must have a bandlimit L, that means that l < L
% For now, it is also one everywhere in < L

function P = kernelP(L, PHI, THETA) 
	addpath('../spharm');


	% The number of sample we have/need:
	M = 0;
	for l = 0:L-1
		M = M + 2*l + 1;
	end
	N = M;
	if nargin < 2
		% - The sampling grid is of size M
		% Define an uniform sampling grid:
		PHI 	=linspace(0, 2*pi - 2*pi/M  , M);  % Azimuthal/Longitude/Circumferential
		THETA  	=linspace(0,   pi - pi/M + 0.001   , M);  % Altitude /Latitude /Elevation
	else
		N = length(THETA);
	end
	

	% Define our matrix P of size MxM
	P = zeros(M,N);


	% For negative polynomial, I'm using this relationship: 
	% http://en.wikipedia.org/wiki/Associated_Legendre_polynomials#Negative_m_and.2For_negative_.E2.84.93
	% So:
	Nneg = @(m, l) (-1)^m * factorial(l-m)/factorial(l+m);

	% For the coefficient, we are using:
	NSpher = @(m, l) sqrt( ((2*l+1)/(4*pi)) *  factorial(l-m)/factorial(l+m)  );


	% Dummy var to know where we are:
	current = 0;
	% fill the alternative matrix with P_l
	for l = 0:L-1
		% Compute the associated legendre polynomials:
		Lmn 	= legendre(l, cos(THETA));
		for i = -l:1:l
			current = current +1;
			if i < 0
				P(current, :) = NSpher( abs(i), l) * Nneg( abs(i), l) .* Lmn( abs(i)+1, :) .* exp( 1i * i * PHI);
			elseif i >= 0
				P(current, :) = NSpher(i, l) .* 		  			    Lmn(      i+1, :) .* exp( 1i * i * PHI);
			end
			% fprintf('constructiong l=%d i=%d, \tY = %g %gi,  \tL = %f, \tNeg = %f\n', l, i, P(current, 1),  P(current, 1),  Lmn( abs(i)+1, 1), NSpher( abs(i), l));
		end


	end

	P = P.';

end

