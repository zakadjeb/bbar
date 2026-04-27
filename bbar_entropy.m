function [H, info] = bbar_entropy(x, m, tau, normalize)
% BBAR_ENTROPY  Permutation entropy of a 1D signal (Bandt & Pompe).
%
%   H = bbar_entropy(x)
%   H = bbar_entropy(x, m, tau, normalize)
%   [H, info] = bbar_entropy(...)
%
% Inputs
% ------
% x          : 1D signal (vector)
% m          : embedding dimension / pattern length (default 5)
% tau        : delay in samples (default 1)
% normalize  : true/false (default true). If true, returns H/log2(m!)
%
% Outputs
% -------
% H          : permutation entropy (bits if normalize=false, else in [0,1])
% info       : struct with fields:
%              .m, .tau, .Npatterns, .counts, .p, .Hraw
%
% Notes
% -----
% - Permutation entropy captures temporal/ordinal structure (rhythm / predictability),
%   unlike histogram-based entropy which ignores time ordering.
% - Ties are handled by adding tiny jitter (stable, deterministic) to break ties.

    if nargin < 2 || isempty(m), m = 5; end
    if nargin < 3 || isempty(tau), tau = 1; end
    if nargin < 4 || isempty(normalize), normalize = true; end

    x = x(:);
    x = x(isfinite(x));
    N = numel(x);

    if N < (m-1)*tau + 1
        error('bbar_entropy:NotEnoughData', ...
            'Signal too short for m=%d, tau=%d (need at least %d samples).', ...
            m, tau, (m-1)*tau + 1);
    end

    % Break ties deterministically (important for EEG where plateaus can occur after rounding)
    % Scale jitter to data range to remain negligible.
    xr = max(x) - min(x);
    if xr == 0
        % constant signal -> only one pattern -> entropy 0
        Hraw = 0;
        H = 0;
        info = struct('m', m, 'tau', tau, 'Npatterns', 1, ...
                      'counts', 1, 'p', 1, 'Hraw', Hraw);
        return;
    end
    jitter = (eps(xr) + 1e-12*xr) * ( (1:N)' - (N+1)/2 ) / N;
    x = x + jitter;

    % Number of ordinal patterns
    nPat = factorial(m);
    counts = zeros(nPat, 1);

    % Precompute Lehmer-code weights to rank permutations in [1..m!]
    % rank = 1 + sum_{i=1..m} c_i * (m-i)! where c_i is number of remaining smaller elements.
    fact = factorial(0:m-1);          % [0!, 1!, ..., (m-1)!]
    fact = fliplr(fact);              % [(m-1)!, ..., 0!]

    % Slide window and count ordinal patterns
    nWindows = N - (m-1)*tau;
    for t = 1:nWindows
        idx = t : tau : (t + (m-1)*tau);
        v = x(idx);

        % Ordinal pattern: permutation of ranks (1..m)
        [~, ord] = sort(v, 'ascend');  % ord is a permutation of 1..m

        % Convert permutation to unique rank using Lehmer code
        % Compute Lehmer code c_i
        c = zeros(1, m);
        remaining = ord;
        for i = 1:m
            c(i) = sum(remaining(1) > remaining(2:end)); %#ok<*BDSCA>
            remaining(1) = [];
        end
        % Rank in 1..m!
        patId = 1 + sum(c .* fact);

        counts(patId) = counts(patId) + 1;
    end

    % Probabilities over observed patterns
    p = counts / sum(counts);
    p = p(p > 0);

    % Shannon entropy of ordinal patterns
    Hraw = -sum(p .* log2(p));  % in bits

    if normalize
        H = Hraw / log2(factorial(m)); % in [0,1]
    else
        H = Hraw;
    end

    info = struct('m', m, 'tau', tau, 'Npatterns', nPat, ...
                  'counts', counts, 'p', p, 'Hraw', Hraw);
end
