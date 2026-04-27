function [stat, p, df, expected, testUsed] = bbar_chisquare(observed, doprint)
% BBAR_CHISQUARE  Chi-square test with automatic fallback to Fisher's exact test
%
%   [stat, p, df, expected, testUsed] = bbar_chisquare(observed, doprint)
%
% Inputs:
%   observed : contingency table
%   doprint  : (optional) 1 = print result (default), 0 = silent
%
% Outputs:
%   stat     : chi-square statistic (or NaN if Fisher used)
%   p        : p-value
%   df       : degrees of freedom (NaN if Fisher used)
%   expected : expected cell counts
%   testUsed : 'chi-square' or 'fisher'

    if nargin < 2
        doprint = 1;
    end

    % Totals
    rowSums = sum(observed, 2);
    colSums = sum(observed, 1);
    grandTotal = sum(observed(:));

    % Expected counts
    expected = (rowSums * colSums) / grandTotal;

    % Check assumptions
    smallExpected = expected < 5;
    tooManySmall = sum(smallExpected(:)) / numel(expected) > 0.2;

    useFisher = any(expected(:) < 1) || tooManySmall;

    % --- CASE 1: Use Fisher (only valid for 2x2) ---
    if useFisher && all(size(observed) == [2 2])
        [~, p] = fishertest(observed);

        stat = NaN;
        df = NaN;
        testUsed = 'fisher';

        if doprint
            fprintf('Fisher''s Exact Test: p = %.4f\n', p);
        end

    % --- CASE 2: Use Chi-square ---
    else
        % Chi-square statistic
        stat = sum((observed - expected).^2 ./ expected, 'all');

        % Degrees of freedom
        df = (size(observed,1)-1) * (size(observed,2)-1);

        % p-value
        p = 1 - chi2cdf(stat, df);

        testUsed = 'chi-square';

        if doprint
            fprintf('Chi-square: χ²(%d) = %.3f, p = %.4f\n', df, stat, p);

            if useFisher && ~all(size(observed) == [2 2])
                fprintf('WARNING: Small expected counts detected, but Fisher test not available for >2x2 tables.\n');
            end
        end
    end
end