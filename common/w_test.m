function w_true = w_test(fri, w)
	% This function compute the true value of w
	% from Eq. (14) of the paper

	K = length(fri.Weights);
	w_true = zeros(K, 1);
	Nn = Normalizing(K+1);

	% Compute the variables:
	c = cos( fri.Locations(:, 1) );
	u = exp(- 1i  .*  fri.Locations(:, 2)  );
	s = sin( fri.Locations(:, 1)  );

	for n = 1:K
	    w_true(n) = sum(  fri.Weights' .*  c .* (s.* u).^(n-1) ) ;
	end


	if nargin > 1
		if abs(w_true - w) < 0.0000000001
			s = 'correct';
		else
			s = 'incorrect';
		end
	  	fprintf('w is %s \n', s);
	end
end