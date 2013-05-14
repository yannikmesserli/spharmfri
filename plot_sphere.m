function  plot_sphere(N, fri)


	figure;

	% generate the sphere:s
	[sx, sy, sz] = sphere(2^4);
	% take only half of it:

	sx = sx(2^(3)+1:end,:);       
	sy = sy(2^(3)+1:end,:);       
	sz = sz(2^(3)+1:end,:);


	surf(sx, sy, sz,  0.5.*ones(size(sz) ) );
	colormap([0 0 0; 0.9 0.9 0.9; 1 1 1]);
	axis equal;
	hold on;


	plotonsphere( fri.Locations(1,1) , fri.Locations(1,2), -3:0.1:3 );
	
	hold off;

end
