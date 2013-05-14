function print_check(name, x, y)

	fprintf(name);

	N = length(x);

	for i = 1:N;
	    fprintf('%-1.4f = %-1.4f', 	x(i), y(i));
	    if i ~= N
	        fprintf(', ');
	    end
	end
	fprintf('\n');

end