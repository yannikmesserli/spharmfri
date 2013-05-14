function s = samples(N, f, plot)
	% Using Eq. (9) of the paper,
	% Plot the Figure 1.
	% kernel is a low-pass filter with bandlimit:
	% L = 2*K-1


	% Sampling grid:
	grid = [500 500];
	% We need only one of the grid:
	Theta = zeros(grid);
	Phi = zeros(grid);

    % The sample:
    s = zeros(grid);


	% First sum, \sum_{l = 0}^L
	for l = 0:N-1
		% Second sum:
        for m = 0:l
            
          	% compute p^_l
            p = zeros(grid);
            % sum of all coeff:
            for n = 0:l
            	% Test if the h= 0 or 1:
                if m*n <= N -1 
	                % Take the spherical harmonics on a grid of theta x phi
	                [Ymn,Theta,Phi,sX,sY,sZ] = spharm(l, m, grid, 0);
	                p  = p + Ymn;
                end
            end
            
            % Compute the sample, Eq. (9):
            s = repmat(f(l+1,m+1), grid) .* p;
            

            
        end
    end

    if plot 

    	
    	figure; title('s(\omega)');
		hold on;
	    % Real part:
	    [X,Y,Z]=sph2cart(Theta,Phi-pi/2,real(s).^2);
	    % Plot it:
	    surf( X, Y, Z, 'EdgeColor','none');


	    axis equal off;
		% %rot3d;
		light; lighting phong; camzoom(1.3);

		hold off;
	end

end


