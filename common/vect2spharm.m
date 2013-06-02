function [m mneg] = vect2spharm(v)
	% Nb of element
	nb = length(v);
	% Find the size of the matrix
	s = sqrt(nb);
	% Test if it's a valid transformation:
		% if ~ isinteger(s)
		% 	throw( MException('vect2spharm:invalide', 'Impossible transformation'));
		% end
	% output:
	m = zeros(s,s);
	mneg = zeros(s,s);
	% Stupid counter
	current = 0;
	for i = 0:s-1
		for j = -i:1:i
			current = current + 1;
			if j == 0
				m(i+1, j+1 ) = v(current);
				mneg(i+1, abs(j)+1 ) = v(current);
			elseif j < 0	
				mneg(i+1, abs(j)+1 ) = v(current);
			else
				m(i+1, j+1 ) = v(current);
			end
		end
	end

end