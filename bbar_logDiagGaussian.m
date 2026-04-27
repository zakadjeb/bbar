function lq = bbar_logDiagGaussian(theta, mu_q, sigma_q)
    dif  = theta - mu_q;
    varq = sigma_q.^2;
    lq   = -0.5 * sum( log(2*pi*varq) + (dif.^2) ./ varq );
end