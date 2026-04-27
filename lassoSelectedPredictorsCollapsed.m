function [selectedDummies, collapsed, intercept] = lassoSelectedPredictorsCollapsed(B, FitInfo, Xnames, tol)
% lassoSelectedPredictorsCollapsed
%
% Inputs:
%   B        : p × nLambda coefficient matrix from lassoglm
%   FitInfo  : FitInfo struct from lassoglm
%   Xnames   : predictor names (string array or cellstr), length p
%   tol      : (optional) non-zero tolerance (default 1e-8)
%
% Outputs:
%   selectedDummies : table of non-zero predictors at lambda_min
%   collapsed       : table collapsing dummy variables into parent factors
%   intercept       : intercept at lambda_min (scalar)
%
% Collapsing rule:
%   Parent factor name = predictor name with the LAST "_<suffix>" removed.
%   Examples: "ASA_1" -> "ASA", "Køn_Mand" -> "Køn", "BMI" -> "BMI"

    if nargin < 4 || isempty(tol)
        tol = 1e-8;
    end

    % Ensure string array
    if iscell(Xnames)
        Xnames = string(Xnames);
    end

    % Basic checks
    if numel(Xnames) ~= size(B,1)
        error("Number of Xnames (%d) must match number of rows in B (%d).", numel(Xnames), size(B,1));
    end
    if ~isfield(FitInfo, "IndexMinDeviance")
        error("FitInfo must contain IndexMinDeviance (use lassoglm with CV).");
    end

    % Pick lambda (min deviance)
    idx = FitInfo.IndexMinDeviance;
    beta = B(:, idx);

    % Intercept at this lambda
    intercept = FitInfo.Intercept(idx);

    % Non-zero (with tolerance)
    keepIdx = abs(beta) > tol;

    % Selected dummy-level predictors
    selectedDummies = table( ...
        Xnames(keepIdx)', ...
        beta(keepIdx), ...
        abs(beta(keepIdx)), ...
        'VariableNames', {'Predictor','Coefficient','AbsCoefficient'} ...
    );

    % If nothing selected, return empty collapsed table too
    if isempty(selectedDummies)
        collapsed = table(string.empty(0,1), zeros(0,1), zeros(0,1), zeros(0,1), strings(0,1), ...
            'VariableNames', {'Factor','NumSelected','MaxAbsCoefficient','SumAbsCoefficient','Members'});
        disp("No predictors survived LASSO at lambda_min (within tolerance).");
        return
    end

    % Compute parent factor names by stripping last "_suffix" if present
    pred = selectedDummies.Predictor;
    parent = regexprep(pred, "_[^_]+$", "");   % remove last underscore segment
    selectedDummies.Factor = parent;

    % Collapse: for each factor, summarize its selected members
    factors = unique(parent, 'stable');

    NumSelected       = zeros(numel(factors),1);
    MaxAbsCoefficient = zeros(numel(factors),1);
    SumAbsCoefficient = zeros(numel(factors),1);
    Members           = strings(numel(factors),1);

    for i = 1:numel(factors)
        f = factors(i);
        m = (parent == f);

        NumSelected(i)       = sum(m);
        MaxAbsCoefficient(i) = max(selectedDummies.AbsCoefficient(m));
        SumAbsCoefficient(i) = sum(selectedDummies.AbsCoefficient(m));

        memNames = selectedDummies.Predictor(m);
        Members(i) = strjoin(memNames, ", ");
    end

    collapsed = table( ...
        factors, ...
        NumSelected, ...
        MaxAbsCoefficient, ...
        SumAbsCoefficient, ...
        Members, ...
        'VariableNames', {'Factor','NumSelected','MaxAbsCoefficient','SumAbsCoefficient','Members'} ...
    );

    % Pretty display (optional)
    disp("Selected predictors (dummy level) at lambda_min:");
    disp(sortrows(selectedDummies, "AbsCoefficient", "descend"));

    disp("Collapsed selection by parent factor:");
    disp(sortrows(collapsed, "MaxAbsCoefficient", "descend"));
end