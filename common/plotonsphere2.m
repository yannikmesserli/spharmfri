function plotonsphere2( theta, phi)
    r = 0:0.1:1;
;
    tot = length(theta);

    hold on;
    for i = 1:tot
        
        x = r .* sin(pi/2 * theta(i)  ) .* cos( 2 * pi * phi(i)  ) ;
        y = r .* sin(pi/2 * theta(i)  ) .* sin( 2 * pi * phi(i)  ) ;
        z = r .* cos(pi/2 * theta(i)  );
        plot3(x, y, z);
        
    end
    hold off;


end

