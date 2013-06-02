function print_check(name, x, y, torerance)

	if nargin < 4
		torerance = 0.000000000001;
	end

	fprintf(name);

	if abs(x - y) < torerance
		s = 'passed! ';
	else
		s = 'failed! ';
	end
	fprintf(s);
	N = length(x);

	for i = 1:N;
	    fprintf('%-1.4f = %-1.4f', 	x(i), y(i));
	    if i ~= N
	        fprintf(', ');
	    end
	end

	% Print the precision:
	%fprintf('(width e=%)')
	
	fprintf('\n');

end