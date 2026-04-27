function effTable = bbar_lmeEffectSize(lme)
% lmeEffectSize  Compute effect sizes for each fixed effect in a LinearMixedModel
%
%   effTable = lmeEffectSize(lme)
%
%   Adds partialEta2 and cohensF for each row of anova(lme).

    % ANOVA table for fixed effects
    a = anova(lme);  % default DFMethod is usually Satterthwaite

    % Extract needed columns
    F   = a.FStat;
    df1 = a.DF1;
    df2 = a.DF2;

    % Preallocate
    partialEta2 = nan(size(F));
    cohensF     = nan(size(F));

    % Compute only where F and dfs are valid (skip intercept if needed)
    valid = ~isnan(F) & df1 > 0 & df2 > 0;

    partialEta2(valid) = (F(valid) .* df1(valid)) ./ ...
                         (F(valid) .* df1(valid) + df2(valid));

    cohensF(valid) = sqrt(partialEta2(valid) ./ (1 - partialEta2(valid)));

    % Add to table
    effTable = a;
    effTable.partialEta2 = partialEta2;
    effTable.cohensF     = cohensF;
end