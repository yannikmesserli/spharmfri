function [h err] = kernelTrain(sample, fri, theta, phi)
% Using the sample and knowing the final result fri, the algorithm will train
% a sampling kernel h so that it minimize the error with fri

	% Parameter:
	K = numel(fri.Weights);
	L = 2*K;
	M = length(sample);
	N = (2*K)^2;


	% Starting value for h:
	h_start = randn(2*(L^2)+1, L+1);
	% function to optimize w.r.t. h
	objFct = @(h) ErrorFor(sample, K, theta, phi, h, fri);
	% Search the min:
	[h err] = fminsearch(objFct, h_start);

end

% Function to optimize:
function err = ErrorFor(sample, K, theta, phi, h, fri)
	% Solve the system
	fri_est = solveFRI(sample, K, theta, phi, h);
	% Compute the error:
	err1 = RMSE(fri_est.Locations, fri.Locations);
	err2 = RMSE(fri_est.Weights, fri.Weights);
	err = mean([ err1 err2 ]);

end