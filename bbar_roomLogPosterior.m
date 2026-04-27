function lp = bbar_LogPosterior(theta, data)
    % theta = [a_raw, log_sigma]
    a_raw     = theta(1);
    log_sigma = theta(2);

    % Transform to constrained parameters
    alpha = 1 ./ (1 + exp(-a_raw));  % (0,1)
    sigma = exp(log_sigma);          % (0,inf)

    % Data
    y  = data.resp;
    HB = data.HB;
    HA = data.HA;

    % Predicted mean: Bayesian-like integration
    mu = alpha .* HB + (1 - alpha) .* HA;

    % Likelihood: y ~ N(mu, sigma^2)
    loglike = sum(log(normpdf(y, mu, sigma)));

    % Priors (same as in MCMC version – weak, regularizing):
    % a_raw ~ N(0, 1)
    % log_sigma ~ N(log(0.1), 1)
    logprior_a    = log(normpdf(a_raw, 0, 1));
    logprior_lsig = log(normpdf(log_sigma, log(0.1), 1));

    lp = loglike + logprior_a + logprior_lsig;
end