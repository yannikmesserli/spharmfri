function r_true = r_test(fri, r)
	% This function compute the true value of w
	% from Eq. (14) of the paper

	K = length(fri.Weights);
	Nn = Normalizing(K+1);

	% Compute the variables:
	u = exp(- 1i  .*  fri.Locations(:, 2)  );
	s = sin( fri.Locations(:, 1)  );

	r_true = u .* s;

	if nargin > 1
		if abs(r_true - r) < 0.0000000001
			s = 'correct';
		else
			s = 'incorrect';
			r
			r_true
		end
	  	fprintf('r is %s \n', s);
	end
end