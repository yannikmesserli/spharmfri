% This function computes the spherical harmonics value 
% at the position theta and phi
% The spherical harmnic order is m and degree is l

function Ylm = harmonic(l, m, theta, phi)
	
	% Check that l > m:
	if l<m, error('The ORDER (m) must be less than or eqaul to the DEGREE(l).'); end


	% Compute the associated legendre polynomials of degree l
	Plm = legendre(l, cos(phi) );
	% Take only the m order:
	Plm = Plm( abs(m), :);
	
	% Handle numerical problem:
	if l~=0
	  Plm = squeeze(Plm(:));
	end

	% Compute the normalizing factor
	a1=((2*l+1)/(4*pi));
	a2=factorial(l - abs(m))/factorial(l+m);
	Nlm = sqrt(a1 * a2);
	
	Ylm =  Plm.' .* exp( i*abs(m)*theta);
	

end