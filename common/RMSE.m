function r = RMSE(y, yhat)
	r = sqrt(mean( abs(y - yhat).^2 ));
end