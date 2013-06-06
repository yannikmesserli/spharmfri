function print_rmse(name_file)

	load(name_file);

	
	plotStyle = {'b','g','r', 'k', 'k:'};
	figure;
	hold on;
	
	[nb_tick blop] = size(rmse_f);

	for i = 1:nb_tick
		% plot(db_axis, rmse_f(i, :), plotStyle{i} );
		% plot(db_axis, rmse_angle(i, :), plotStyle{i} );
		plot(db_axis, rmse_angle(i, :), plotStyle{i} );
	end

	legend(legends);
	xlabel('SNR in dB');
	ylabel('RMSE');

	hold off;

end