function [x_fit, y_fit, mdl, offset, cubic] = bbar_cubic_fitting(x,y)
% Fits y = b0 + b3*x^3 
%
% Returns:
%   x_fit, y_fit : smooth fitted curve over data range
%   mdl          : fitnlm model (use coefCI, predict, etc.)
%   offset, linear, cubic : estimated coefficients b0, b1, b3

    % --- tidy inputs ---
    x = x(:); y = y(:);
    ok = isfinite(x) & isfinite(y);
    x = x(ok); y = y(ok);

    % --- model ---
    % cubicModel = @(b,x) b(1) + b(2).*x.^2 + b(3).*x.^3;  % [b0 b1 b3]
    cubicModel = @(b,x) b(1) + b(2).*x.^3; % + b(3).*x.^3;  % [b0 b1 b3]

    % --- initial guess via quick linear least squares on [1, x, x^3] ---
    X0 = [ones(size(x)) x.^3];
    beta0_ls = X0 \ y;
    if any(~isfinite(beta0_ls))
        beta0_ls = [median(y); 0; 0]; % fallback
    end
    beta0 = beta0_ls.';  % row vector as fitnlm expects

    % --- fit with fitnlm ---
    tbl = table(x,y);
    mdl = fitnlm(tbl, cubicModel, beta0);

    % --- output params ---
    p = mdl.Coefficients.Estimate;
    offset = p(1);
    % linear = p(2);
    cubic  = p(2);

    % --- smooth curve for plotting ---
    x_fit = linspace(min(x), max(x), 500).';
    y_fit = cubicModel(p, x_fit);
end
