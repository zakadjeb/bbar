function out = bbar_sigmoid_fitting(x, y, group, V)
% NLME Gompertz with offset:
%   y = d + (a - d) * exp( - b * exp( - c * x ) )
%
% Inputs:
%   x      Nx1 predictor
%   y      Nx1 response
%   group  Nx1 grouping variable (numeric, string/char, or categorical)
%   V      MxG group-level predictors (or []), M = number of groups in grp2idx(group)
%
% Output struct 'out':
%   beta, PSI, stats, B   : nlmefitsa outputs (fixed, RE cov, stats, RE per group)
%   s                     : scale used on x (x_s = x/s)
%   phi_pop               : population parameters [a b c d] (on original scale)
%   phi_subject(i,:)      : subject-specific params in grp2idx order
%   predict(xnew, subj)   : prediction handle (subj=[] → population)
%   ciFixed               : 95% Wald CI for fixed effects on *transformed* scale

    % --- tidy ---
    x = x(:); y = y(:);
    ok = isfinite(x) & isfinite(y);
    x = x(ok); y = y(ok);
    group = group(ok);

    % --- group mapping (defines row order for V, B, etc.) ---
    [gid, grpNames] = grp2idx(group);
    M = max(gid);

    % --- scale x (no centering) for numerics ---
    s  = max(abs(x)); if s==0, s=1; end
    xs = x / s;
    X  = xs;                  % predictor matrix for nlmefitsa (N×H). Here H=1.

    % --- V: group-level predictors (may be []) ---
    if nargin < 4 || isempty(V)
        Vuse = [];
    else
        assert(size(V,1)==M, 'V must have M rows in grp2idx(group) order.');
        Vuse = V;
    end

    % --- model function: MODELFUN(PHI, Xfun, Vfun) -> yhat (vector)
    % PHI is 1×4 here: [a b c d]
    modelfun = @(phi, Xfun, Vfun) ...
        phi(4) + (phi(1) - phi(4)) .* exp( - phi(2) .* exp( - phi(3) .* Xfun(:,1) ) );

    % --- starting values (on original scale) ---
    d0 = min(y);
    a0 = max(y) - d0;
    c0 = 1 / max(std(xs), eps);  % growth rate on scaled x
    b0 = 1;

    % We'll log-transform a,b,c to enforce positivity (d untransformed).
    % With ParamTransform=[1 1 1 0], BETA are fixed effects on the *log* scale for a,b,c.
    beta0_orig = [a0, b0, c0, d0];     % [a b c d] on original scale
    beta0      = [log(max(a0,eps)); log(max(b0,eps)); log(max(c0,eps)); d0];

    xform = [1 1 1 0];  % 1=log, 0=identity (for [a b c d])

    % Random effects: by default, one per parameter (diagonal PSI). You can adjust via options.
    reSel = [true false false true];

    % --- fit ---
    [beta, PSI, stats, B] = nlmefitsa( ...
        X, y, gid, Vuse, modelfun, beta0, ...
        'ParamTransform', xform, ...
        'REParamsSelect', reSel, ...
        'CovPattern', eye(sum(reSel)), ...
        'ErrorModel', 'constant', ...
        'ComputeStdErrors', true);

    % --- helper: inverse transform to get actual PHI from fixed (and random) effects
    % On transformed scale: theta = beta + b_i (default FE/RE design = I)
    invT = @(theta) [exp(theta(1)), exp(theta(2)), exp(theta(3)), theta(4)];

    % Population parameters (random = 0)
    phi_pop = invT(beta(:).');         % 1×4 [a b c d] on original scale

    % Subject-specific parameters in grp2idx order
    phi_subject = zeros(M,4);
    for i = 1:M
        phi_subject(i,:) = invT(beta(:).' + B(:,i).');
    end

    % --- pack outputs ---
    out.beta        = beta;
    out.PSI         = PSI;
    out.stats       = stats;
    out.B           = B;
    out.s           = s;
    out.grpNames    = grpNames;
    out.phi_pop     = phi_pop;
    out.phi_subject = phi_subject;

    % prediction handles (xnew on original x-scale)
    out.predict = @(xnew, subj) local_predict(xnew, subj, s, modelfun, phi_pop, grpNames, phi_subject);
    out.ciFixed = @() wald_ci_fixed(stats, beta);   % CI for beta on transformed scale
end

% ---- prediction (population or subject-specific)
function yhat = local_predict(xnew, subj, s, modelfun, phi_pop, grpNames, phi_subject)
    xsnew = xnew(:) / s;
    if isempty(subj)
        phi = phi_pop;
    else
        % subj can be index in grp2idx order, or a label matching grpNames
        if isnumeric(subj) && isscalar(subj) && subj>=1 && subj<=size(phi_subject,1)
            phi = phi_subject(subj,:);
        else
            % match by name
            if ischar(subj) || isstring(subj), subj = string(subj); end
            idx = find(strcmp(grpNames, cellstr(subj)), 1);
            if isempty(idx), error('Subject not found in grp2idx(group) order.'); end
            phi = phi_subject(idx,:);
        end
    end
    yhat = modelfun(phi, xsnew, []);
end

% ---- Wald 95% CI for fixed effects (on transformed scale!)
function CI = wald_ci_fixed(stats, beta)
    if isfield(stats,'covb') && ~isempty(stats.covb)
        se = sqrt(diag(stats.covb));
        z  = 1.95996398454005;
        CI = [beta(:)-z*se, beta(:)+z*se];
    else
        CI = NaN(numel(beta),2);
    end
end
