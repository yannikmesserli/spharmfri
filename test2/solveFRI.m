function fri_est = solveFRI(samples, K, phi, theta)

	
	

	nb_samples = length(samples);

	N = 2*K;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% COMPUTE THE SPHERICAL HARMONICS FROM THE SAMPLES:
	%

	% Construct the matrix P
	% Which consist of only the independant column of the 
	% kernel of a low-pass filter:
	

	if nargin < 3
		% The number of sample we have/need:
	else
		P = kernelP(K);
	end

	% The coefficient of same order and same degree
	f_est = pinv(P) * samples;
	% Put it back on my standart form:
	[f fNeg] = vect2spharm(f_est);
	% Catch only the diagonal:
	fnn = diag(f).';
	% Catch only the diagonal below
	fn1n = diag(f, -1)';

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% FIND PHI:
	
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
	

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% FIND THE AMPLITUDE:
	%

	% Create the Vandermonde matrix:
	V = repmat( r, 1, K).'  .^ repmat( (0:K-1).', 1, K);

	% Solve it:
	a = sort( real(V \ z(1:K).'));

	

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

	

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% CONSTRUCT THE ESTIMATED SIGNAL:
	%
	
	fri_est.Locations = [  theta_est phi_est];
	fri_est.Weights = a';

end