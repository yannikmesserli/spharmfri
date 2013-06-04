function plotonsphere( theta, phi, r )
        
        hold on;
        x = r .* sin(pi/2 * theta  ) .* cos( 2 * pi * phi  ) ;
        y = r .* sin(pi/2 * theta  ) .* sin( 2 * pi * phi  ) ;
        z = r .* cos(pi/2 * theta  );
        plot3(x, y, z);
        hold off;


end

