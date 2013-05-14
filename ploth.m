% Plot the sampling kernel h
addpath('spharm');

figure; title('Sampling kernel:');

X = zeros(100, 100);
Y = zeros(100, 100);
Z = zeros(100, 100);
for l = 0:4
    for n = 0:l
        [Ymn,THETA,PHI,sX,sY,sZ] = spharm(l, n, [100 100], 0);



         % abs part:

        [Xr,Yr,Zr]=sph2cart(THETA,PHI-pi/2,real(Ymn).^2);
        X = X + Xr;
        Y = Y + Yr;
        Z = Z + Zr;
        
    end
end

hold on;
surf( X, Y, Z);
%axis equal off; %rot3d;
hold off;
light; lighting phong; camzoom(1.3);