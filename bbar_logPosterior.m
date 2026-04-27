function lp = bbar_logPosterior(theta, data)
    % theta = [beta0, beta1, log_sigma0, log_gamma]
    beta0      = theta(1);
    beta1      = theta(2);
    log_sigma0 = theta(3);
    log_gamma  = theta(4);

    sigma0 = max(exp(log_sigma0), 1e-6);
    gamma  = max(exp(log_gamma),  1e-6);

    y     = data.y;
    HB    = data.HB;
    ratio = data.ratio;

    % Bias and mean per trial
    bias_i = beta0 + beta1 .* (ratio - 1);
    mu_i   = HB + bias_i;

    % Trial-wise variance
    var_i = sigma0.^2 + gamma .* (ratio - 1).^2;
    sd_i  = sqrt(var_i);

    % Likelihood: y_i ~ N(mu_i, var_i)
    loglike = sum(log(normpdf(y, mu_i, sd_i)));

    % Priors:
    % beta0 ~ N(0, 1)
    % beta1 ~ N(0, 1)
    % log_sigma0 ~ N(log(0.1), 1)
    % log_gamma  ~ N(log(0.1), 1)
    logprior_b0     = log(normpdf(beta0, 0, 1));
    logprior_b1     = log(normpdf(beta1, 0, 1));
    logprior_sigma0 = log(normpdf(log_sigma0, log(0.1), 1));
    logprior_gamma  = log(normpdf(log_gamma,  log(0.1), 1));

    lp = loglike + logprior_b0 + logprior_b1 + logprior_sigma0 + logprior_gamma;
end