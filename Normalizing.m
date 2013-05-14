function Nn = Normalizing(N)

	% Computed experimentally:
	%NBig = [1.0000   -1.0000    3.0000  -15.0000  105.0000 -945.0000];
	% Nn = NBig(1:N);

	n = 1:N-1;
	Nn = ones(N,1);
	Nn(2:end) = (-1).^n .* factorial(2*n-2) ./ factorial(n-1) .* (n - 0.5)./2.^(n-2);

	% Comment:
	% WTF?????:
end