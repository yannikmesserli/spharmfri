% We are computing the the matrix P using Eq. (8) in the paper. 
% - Our kernel must have a bandlimit L, that means that l < L
% For now, it is also one everywhere in < L

function [P Palt] = kernelP(L) 
	addpath('../spharm');
	
	% Define our matrix P of size MxM
	P = zeros(L,L);

	
	% - The sampling grid is of size M
	% Define an uniform sampling grid:
	THETA	=linspace(0, 2*pi, L);  % Azimuthal/Longitude/Circumferential
	PHI  	=linspace(0,   pi, L);  % Altitude /Latitude /Elevation

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
			P(l+1, : ) =  P(l+1, : ) + NSpher(0, l) .* Lmn( 1, :);
		end

		% Compute the sum Eq. (8)
		for i = -l:1:l
			if i < 0
				% Compute for the negative one, add the coeff:
				P(l+1, : ) =  P(l+1, : ) + NSpher(i, l) * Neq(i, l) .* Lmn( abs(i)+1, :) .* exp(1i*i*THETA);
			else
				P(l+1, : ) =  P(l+1, : ) + NSpher(i, l) .* 			 Lmn(      i+1, :) .* exp(1i*i*THETA);
			end


			
		end



		
	end

end




