function [smd,labels] = compute_smd(x, treat, w)
% COMPUTE_SMD Computes standardized mean difference
% Handles numeric, logical, binary categorical, and multi-level categorical variables.
% Ignores rows with invalid weights (NaN, Inf, negative), missing treatment,
% and missing x values where relevant.
%
% Inputs:
%   x     : covariate
%   treat : 0/1 treatment indicator
%   w     : weights
%
% Output:
%   smd   : scalar SMD for numeric/binary variables
%           vector of SMDs for multi-level categorical variables

    if nargin < 3 || isempty(w)
        w = ones(size(treat));
    end

    % Ensure column vectors
    x = x(:);
    treat = treat(:);
    w = w(:);

    % Basic validity for treatment and weights
    valid = ~isnan(treat) & ~isnan(w) & isfinite(w) & (w >= 0);

    x = x(valid);
    treat = treat(valid);
    w = w(valid);

    % Split treatment groups
    idx_t = (treat == 1);
    idx_c = (treat == 0);

    % ---- CASE 1: numeric or logical ----
    if isnumeric(x) || islogical(x)

        x = double(x);

        % Remove NaNs in x separately
        valid_x = ~isnan(x);
        x = x(valid_x);
        treat_num = treat(valid_x);
        w_num = w(valid_x);

        idx_t = (treat_num == 1);
        idx_c = (treat_num == 0);

        x_t = x(idx_t);
        x_c = x(idx_c);
        w_t = w_num(idx_t);
        w_c = w_num(idx_c);
        labels = categorical("Continuous");

        % Guard against empty groups or zero total weight
        if isempty(x_t) || isempty(x_c) || sum(w_t) == 0 || sum(w_c) == 0
            smd = NaN;
            return
        end

        % Weighted means
        m_t = sum(w_t .* x_t) / sum(w_t);
        m_c = sum(w_c .* x_c) / sum(w_c);

        % Weighted variances
        v_t = sum(w_t .* (x_t - m_t).^2) / sum(w_t);
        v_c = sum(w_c .* (x_c - m_c).^2) / sum(w_c);

        denom = sqrt((v_t + v_c) / 2);

        if denom == 0 || isnan(denom)
            smd = 0;
        else
            smd = (m_t - m_c) / denom;
        end

        return
    end

    % ---- CASE 2: categorical / string / char ----
    if iscategorical(x) || isstring(x) || iscellstr(x) || ischar(x)

        % Convert to categorical for consistent handling
        x = categorical(x);

        % Remove missing categories
        valid_x = ~ismissing(x);
        x = x(valid_x);
        labels = unique(x);
        treat_cat = treat(valid_x);
        w_cat = w(valid_x);

        idx_t = (treat_cat == 1);
        idx_c = (treat_cat == 0);

        levels = categories(x);
        smd = NaN(numel(levels), 1);

        for k = 1:numel(levels)
            x_bin = double(x == levels{k});

            x_t = x_bin(idx_t);
            x_c = x_bin(idx_c);
            w_t = w_cat(idx_t);
            w_c = w_cat(idx_c);

            % Guard against empty groups or zero total weight
            if isempty(x_t) || isempty(x_c) || sum(w_t) == 0 || sum(w_c) == 0
                smd(k) = NaN;
                continue
            end

            % Weighted proportions
            p_t = sum(w_t .* x_t) / sum(w_t);
            p_c = sum(w_c .* x_c) / sum(w_c);

            % Pooled variance for binary variable
            v_t = p_t * (1 - p_t);
            v_c = p_c * (1 - p_c);

            denom = sqrt((v_t + v_c) / 2);

            if denom == 0 || isnan(denom)
                smd(k) = 0;
            else
                smd(k) = (p_t - p_c) / denom;
            end
        end

        return
    end

    error('Unsupported variable type for SMD computation.');
end