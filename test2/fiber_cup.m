% In this file, we are extracting the samples from the fiber cup data
% http://www.lnao.fr/spip.php?rubrique79

clear all;


% We will 
addpath('../common');
addpath('../nifti');

% Read the file containing the direction
directions = dlmread('../fiber_cup/diffusion_directions_corrected.txt');
% Get the sampled position:
[phi theta r] = cart2sph(directions);

figure; plotonsphere2(phi, theta);

% power = 2;
% s_noisy = s + 10^(-power) .* randn(size(s)) + 10^(-power) .* randn(size(s)) * 1i;
% s_noisy

% f_noisy =  pinv(P) * s_noisy;

% abs_diff_f = mean( abs(  (f_noisy - f_true)   ) )
% abs_diff_s = mean( abs(  (s_noisy - s)   ) )

% Then solve it:
fri_est = solveFRI(s, 2);
