function fri_est = solveFRI(samples, fri)

	check = 1;
	if nargin < 2
		% The number of sample we have/need:
		check = 0;
	end

	nb_samples = length(samples);

	K = (nb_samples + 1)/4;
	N = 2*K;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% COMPUTE THE SPHERICAL HARMONICS FROM THE SAMPLES:
	%

	% Construct the matrix P
	% Which consist of only the independant column of the 
	% kernel of a low-pass filter:
	P = kernelP(N);

	% The coefficient of same order and same degree
	fnn = pinv(P) * samples(1:N)';
	fnn = fnn.';
	% The coefficient of degree one minus the order
	fn1n = pinv(P(2:end, 2:end)) * samples(N+1:end)';
	fn1n = fn1n.';

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% FIND PHI:
	%

	Nn = Normalizing(N);
	% Construct z:
	z = fnn ./ Nn;


	% Annihilating filter:
	X = toeplitz( z(K:N-1), fliplr(  z(1:K) )  );
	Y = - z(K+1:N).';
	A = X \ Y;
	r = roots([1 A.']);

	% We can find phi:
	[phi_est, IX] = sort( mod( - angle(r), 2*pi));


	% Put the r_k vector in the right order:
	r = r(IX);
	% Check if we have the correct order:
	if check
		r_test(fri, (r));
		print_check('Position check on phi: ',  phi_est, fri.Locations(:, 2));
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% FIND THE AMPLITUDE:
	%

	% Create the Vandermonde matrix:
	V = repmat( r, 1, K).'  .^ repmat( (0:K-1).', 1, K);

	% Solve it:
	a = sort( real(V \ z(1:K).'));

	if check
		print_check('Amplitude check: ',  a', fri.Weights);
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% FIND THETA:
	% Annihilating filter:
	%

	% Compute w from Eq. (14):
	w = -1 .* fn1n(1:K) ./  Nn(2:K+1);
	w = w.';

	% Check if the vector w is correct:
	if check
		w_test(fri, w);
	end


	% Create the Vandermonde matrix:
	W = repmat( r.', K, 1)  .^ repmat( (0:K-1)', 1, K) .* repmat( a', K, 1);


	% Solve it:
	c_est = (W \ w);

	% We can find theta:
	theta_est   =  sort(mod( real(acos(c_est)), pi));

	if check
		print_check('Position check on theta: ', theta_est, fri.Locations(:, 1));
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% CONSTRUCT THE ESTIMATED SIGNAL:
	%
	
	fri_est.Locations = [  theta_est phi_est];
	fri_est.Weights = a';

end