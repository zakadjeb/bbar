function nELBO = bbar_negELBO(phi, data, eps_mat, D)
    % phi = [mu_q(1:D), log_sigma_q(1:D)]
    mu_q      = phi(1:D);
    log_sig_q = phi(D+1:end);
    sigma_q   = exp(log_sig_q);

    % Reparameterization with fixed eps_mat
    theta = mu_q' + (sigma_q' .* eps_mat');  % D x nSamples
    theta = theta';                          % nSamples x D

    nSamples = size(theta,1);
    logw = zeros(nSamples,1);

    for k = 1:nSamples
        th      = theta(k,:);                             % [beta0, beta1, log_sigma0, log_gamma]
        log_p   = bbar_logPosterior(th, data);       % log p(theta, data)
        log_q   = bbar_logDiagGaussian(th, mu_q, sigma_q);     % log q(theta)
        logw(k) = log_p - log_q;
    end

    ELBO  = mean(logw);
    nELBO = -ELBO;
end