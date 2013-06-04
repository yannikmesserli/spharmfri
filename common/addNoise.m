function s = addNoise(s0, SNRdb)

	P0 = mean(s0.^2);

	SNR = 10^(SNRdb/10);

	s = s0 + sqrt(P0/SNR)*randn(size(s0));

end