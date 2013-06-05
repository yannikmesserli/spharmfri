function print_rmse(name_file)

	load(name_file);

	
	plotStyle = {'b','g','r', 'k', 'k:'};
	figure;
	hold on;
	
	
	%plot(db_axis, rmse_l(i, :), plotStyle{i} );
	% plot(db_axis, rmse_a(i, :), plotStyle{i} );
	plot(db_axis, rmse_f, plotStyle{1} );



	xlabel('SNR in dB');
	ylabel('RMSE');

	hold off;

end