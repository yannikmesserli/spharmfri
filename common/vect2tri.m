function m = vect2tri(v)

	nb = length(v);
	% size of the matrix:
	s = floor(sqrt(2*nb))
	% Test if it's a valid operation
	%if isinteger( s*(s+1)/2 )
		% output:
		m = zeros(s, s);
		% Stupid counter
		current = 0;
		for i = 1:s
			for j = 1:i
				current = current + 1;
				 m(i, j) = v(current);
			end
		end
	% else
	% 	% Error, impossible to do transform it 
	% 	throw( MException('vect2tri:invalide', 'Impossible transformation'));
	% end

end