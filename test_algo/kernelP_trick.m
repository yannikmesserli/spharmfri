% We are computing the the matrix P using Eq. (8) in the paper. 
% - Our kernel must have a bandlimit L, that means that l < L
% For now, it is also one everywhere in < L

function [P Palt] = kernelP(L, M, h) 
	addpath('./spharm');
	% later define h
	% we could use h to search min...


	if nargin < 2
		% The number of sample we have/need:
		M = 0;
		for l = 0:L-1
			M = M + 2*l + 1;
		end
	end

	% Define our matrix P of size MxM
	P = zeros(M,M);

	% Define an alternative matrix P that does not depends on m
	Palt = zeros(L,M);

	% - The sampling grid is of size M
	% Define an uniform sampling grid:
	THETA	=linspace(0, 2*pi, M);  % Azimuthal/Longitude/Circumferential
	PHI  	=linspace(0,   pi, M);  % Altitude /Latitude /Elevation

	% For negative polynomial, I'm using this relationship: 
	% http://en.wikipedia.org/wiki/Associated_Legendre_polynomials#Negative_m_and.2For_negative_.E2.84.93
	% So:
	Neq = @(m, l) (-1)^m * factorial(l-m)/factorial(l+m);

	% For the coefficient, we are using:
	NSpher = @(m, l) sqrt( ((2*l+1)/(4*pi)) *  factorial(l-m)/factorial(l+m)  );


	% Dummy var to know where we are:
	current = 0;

	% fill the alternative matrix with P_l
	for l = 0:L-1

		% Compute the associated legendre polynomials:
		Lmn=legendre(l, cos(PHI));
		%Lmn=squeeze(Lmn(M+1,:,:));

		if l == 0
			Palt(l+1, : ) =  Palt(l+1, : ) + NSpher(0, l) .* Lmn( 1, :);
		end

		% Compute the sum Eq. (8)
		for i = -l:1:l
			if i < 0
				% Compute for the negative one, add the coeff:
				Palt(l+1, : ) =  Palt(l+1, : ) + NSpher(i, l) * Neq(i, l) .* Lmn( abs(i)+1, :) .* exp(1i*i*THETA);
			else
				Palt(l+1, : ) =  Palt(l+1, : ) + NSpher(i, l) * 			 Lmn(      i+1, :) .* exp(1i*i*THETA);
			end


			
		end

		% Fill the P matrix by replicating the value:
		% Update where we have to fill:
		current = current + max( 2*(l - 1)+1, 0);
		nbrep = 2*l + 1;

		P(current+1:current + nbrep, :) = repmat(Palt(l+1, :), nbrep, 1 );


		
	end

end




