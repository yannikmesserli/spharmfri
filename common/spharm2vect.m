function v = spharm2vect(m, mneq)
	% Nb of element
	s = min(size(m));
	% Test that the second matrix is of same size:
	if min(size(mneq)) ~= s
		throw( MException('spharm2vect:invalide', 'Impossible transformation'));
	end
	% Number of element of the lower triangle part of m and the lower triangle part minus the diagonal of mneq
	nb = s*(s+1)/2 + s*(s-1)/2; % = s*s; 
	% output:
	v = zeros(nb, 1);
	% Stupid counter
	current = 0;
	for i = 0:s-1
		for j = -i:1:i
			current = current + 1;
			if j < 0
				v(current) = mneq(i+1, abs(j)+1 );
			else
				v(current) = m(i+1, j+1 );
			end
		end
	end

end