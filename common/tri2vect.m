function v = tri2vect(m)

	s = min(size(m));
	% Number of element of the lower triangle part of m
	nb = s*(s+1)/2;
	% output:
	v = zeros(nb, 1);
	% Stupid counter
	current = 0;
	for i = 1:s
		for j = 1:i
			current = current + 1;
			v(current) = m(i, j);
		end
	end

end