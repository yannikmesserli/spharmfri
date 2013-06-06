close all;
clear all;

k = 4; % definition of the sphere
nb = 5; % nb of samples on theta or phi, better if odd...

% generate the sphere:s
[sx, sy, sz] = sphere(2^k);
% take only half of it:

sx = sx(2^(k-1)+1:end,:);       
sy = sy(2^(k-1)+1:end,:);       
sz = sz(2^(k-1)+1:end,:); 

% size of the samples:


% plot the sphere:
surf(sx, sy, sz,  0.5.*ones(size(sz) ) );
colormap([0 0 0; 0.9 0.9 0.9; 1 1 1]);
axis equal;
hold on;


vect = 0:1/10:1;
for phi = 0:(nb*2)
    for theta = 0:ceil(nb / 2)
        r = (rand(1) * vect) + 1;
        x = r .* sin(pi/2 * theta / ceil(nb / 2) ) .* cos( 2 * pi * phi / (nb *2) ) ;
        y = r .* sin(pi/2 * theta / ceil(nb / 2) ) .* sin( 2 * pi * phi / (nb *2) ) ;
        z = r .* cos(pi/2 * theta / ceil(nb / 2) );
        plot3(x, y, z);
    end
end


