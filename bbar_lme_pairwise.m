function [results, cmap] = bbar_lme_pairwise(model, predictorName, varargin)
% BBAR_LME_PAIRWISE Computes pairwise comparisons and Cohen's d for a categorical
% predictor in a fitlme model, with raincloud-style plotting.
% Also handles categorical:continuous interactions using Johnson-Neyman technique.
%
% INPUTS:
%   model         - fitted fitlme model object
%   predictorName - string, name of the categorical predictor
%
% OPTIONAL PARAMETERS (Name-Value pairs):
%   'alpha'           - transparency level for scatter and density plots (default: 0.5)
%   'plot'            - whether to generate a plot (default: true)
%   'boxplot'         - show boxplot (default: 1)
%   'scatter'         - show scatter plot (default: 1)
%   'normal'          - show kernel density curve (default: 1)
%   'BoxWidth'        - width of boxplot (default: 0.25)
%   'BoxFaceAlpha'    - boxplot face transparency (default: 0.5)
%   'jitterAmount'    - horizontal jitter amount for scatter points (default: 0.1)
%   'ScatterSize'     - size of scatter points (default: 15)
%   'Bandwidth'       - bandwidth for kernel density estimation (default: [], auto)
%   'showMean'        - overlay mean line across conditions (default: 1)
%   'showMedian'      - overlay median line across conditions (default: 1)
%   'ylabel'          - y-axis label (default: 'Magnitude')
%   'xlabel'          - x-axis label (default: '')
%   'xtickslabels'    - cell array of custom x-tick labels (default: {}, uses model level names)
%   'title'           - plot title (default: 'Pairwise comparison')
%   'fontName'        - font name for all text (default: 'Arial')
%   'fontSize'        - font size for all text (default: 10)
%   'grid'            - grid visibility: 'on' or 'off' (default: 'on')
%   'bgColor'         - figure background color (default: 'w')
%   'printP'          - print exact p-values alongside significance stars (default: 1)
%   'StatsYPos'       - [ypos, spacing] manually set position of significance brackets (default: [])
%   'fPosition'       - figure position and size [x y width height] (default: [200 200 600 400])
%   'nexttile'        - plot into existing tiledlayout tile; 0 = new figure (default: 0)
%   'doPlotJN'        - plot Johnson-Neyman regions (default: 0)
%   'nPoints'         - number of points along continuous variable for JN plot (default: 200)
%   'sigAlpha'        - significance level for Johnson-Neyman (default: 0.05)
%
% OUTPUTS:
%   results - table with pairwise comparison statistics
%   cmap    - colormap used in plot
%
% Zakaria Djebbara, March 2026

% --- Detect interaction syntax ---
isInteraction = contains(predictorName, ':');

if isInteraction
    parts = strsplit(predictorName, ':');
    predA = parts{1};
    predB = parts{2};
else
    predA = predictorName;
    predB = '';
end

% --- Input Parser ---
p = inputParser;
addRequired(p, 'model');
addRequired(p, 'predictorName', @ischar);
addParameter(p, 'alpha', 0.5, @isnumeric);
addParameter(p, 'plot', true, @islogical);
addParameter(p, 'boxplot', 1, @isnumeric);
addParameter(p, 'scatter', 1, @isnumeric);
addParameter(p, 'normal', 1, @isnumeric);
addParameter(p, 'BoxWidth', 0.25, @isnumeric);
addParameter(p, 'BoxFaceAlpha', 0.5, @isnumeric);
addParameter(p, 'jitterAmount', 0.1, @isnumeric);
addParameter(p, 'showMean', 1, @isnumeric);
addParameter(p, 'showMedian', 1, @isnumeric);
addParameter(p, 'ylabel', 'Magnitude', @ischar);
addParameter(p, 'xlabel', '', @ischar);
addParameter(p, 'xtickslabels', {}, @iscell);
addParameter(p, 'title', ['Pairwise comparison'], @ischar);
addParameter(p, 'fontName', 'Arial', @ischar);
addParameter(p, 'fontSize', 10, @isnumeric);
addParameter(p, 'grid', 'on', @ischar);
addParameter(p, 'bgColor', 'w');
addParameter(p, 'printP', 1, @isnumeric);
addParameter(p, 'StatsYPos', [], @isnumeric);
addParameter(p, 'fPosition', [200 200 600 400], @isvector);
addParameter(p, 'nexttile', 0);
addParameter(p, 'ScatterSize', 15, @isnumeric);
addParameter(p, 'Bandwidth', [], @isnumeric);
addParameter(p, 'nPoints', 200, @isnumeric);
addParameter(p, 'sigAlpha', 0.05, @isnumeric);
addParameter(p, 'doPlotJN', 0, @isnumeric);
addParameter(p, 'cmap', [], @isnumeric);

parse(p, model, predictorName, varargin{:});

doPlot      = p.Results.plot;

% --- Extract model info ---
coeffNames = model.CoefficientNames';
fe         = fixedEffects(model);
sigma      = sqrt(model.MSE);

% --- Detect whether interaction is cat x cont or simple --- 
if isInteraction
    varA = model.Variables.(predA);
    varB = model.Variables.(predB);
    predA_isCat = iscategorical(varA) || ischar(varA) || iscell(varA) || islogical(varA);
    predB_isCat = iscategorical(varB) || ischar(varB) || iscell(varB) || islogical(varB);
    isCatCont   = predA_isCat && ~predB_isCat;
    isContCat   = ~predA_isCat && predB_isCat;

    % Ensure predA is always the categorical one
    if isContCat
        tmp   = predA; predA = predB; predB = tmp;
        isCatCont = true;
    end

    if ~isCatCont
        error('bbar_lme_pairwise currently supports categorical:continuous interactions only for the JN technique.');
    end
end

% --- Colormap ---
levelsA = categories(categorical(model.Variables.(predA)));
nA      = length(levelsA);

if ~isempty(p.Results.cmap)
    cmap = p.Results.cmap;
else
    cmap = bbar_makecmap([51/255, 92/255, 103/255; 224/255, 159/255, 62/255; 84/255, 11/255, 14/255], nA);
end

% CASE 1: Simple predictor (no interaction)
% =========================================================
if ~isInteraction

    idx    = find(contains(coeffNames, predA) & ~contains(coeffNames, ':'));
    nPairs = nchoosek(nA, 2);

    Group1     = cell(nPairs, 1);
    Group2     = cell(nPairs, 1);
    MeanDiff   = zeros(nPairs, 1);
    SE         = zeros(nPairs, 1);
    tStat      = zeros(nPairs, 1);
    DF         = zeros(nPairs, 1);
    PValue     = zeros(nPairs, 1);
    CohensD    = zeros(nPairs, 1);
    compMatrix = zeros(nPairs, 2);

    row = 1;
    for i = 1:nA
        for j = i+1:nA
            C = zeros(1, model.NumCoefficients);
            if i > 1; C(idx(i-1)) =  1; end
            if j > 1; C(idx(j-1)) = -1; end

            [p_val, F, ~, df2] = coefTest(model, C);
            t = sqrt(F) * sign(C * fe);

            Group1{row}       = levelsA{i};
            Group2{row}       = levelsA{j};
            MeanDiff(row)     = C * fe;
            SE(row)           = sqrt(C * model.CoefficientCovariance * C');
            tStat(row)        = t;
            DF(row)           = df2;
            PValue(row)       = p_val;
            CohensD(row)      = MeanDiff(row) / sigma;
            compMatrix(row,:) = [i, j];
            row = row + 1;
        end
    end

    results = table(Group1, Group2, MeanDiff, SE, tStat, DF, PValue, CohensD);
    fprintf('\nPairwise comparisons for: %s\n', predA);
    fprintf('Effect size (Cohen''s d) based on model residual SD\n\n');
    disp(results);

    if ~doPlot; return; end

    rawData = cell(1, nA);
    for i = 1:nA
        mask = strcmp(string(model.Variables.(predA)), levelsA{i});
        rawData{i} = model.Variables.(model.ResponseName)(mask);
    end

    if p.Results.nexttile == 0
        figure('Position', p.Results.fPosition);
    elseif strcmp(p.Results.nexttile, 'on')
        nexttile;
    else
        nexttile(p.Results.nexttile);
    end

    hold on;
    jitter = bbar_remapnumbers(rand(1, max(cellfun(@length, rawData))), ...
        -p.Results.jitterAmount, p.Results.jitterAmount);

    if p.Results.boxplot
        for i = 1:nA
            boxchart(repmat(i - 0.1, length(rawData{i}), 1), rawData{i}, ...
                'MarkerStyle', 'none', 'BoxFaceColor', cmap(i,:), ...
                'BoxWidth', p.Results.BoxWidth, 'BoxFaceAlpha', p.Results.BoxFaceAlpha);
        end
    end
    if p.Results.scatter
        for i = 1:nA
            n = length(rawData{i});
            scatter(i + jitter(1:n), rawData{i}, p.Results.ScatterSize, cmap(i,:), ...
                'filled', 'MarkerFaceAlpha', p.Results.alpha);
        end
    end
    if p.Results.normal
        for i = 1:nA
            [f, xi] = ksdensity(rawData{i}, 'Bandwidth', p.Results.Bandwidth);
            patch(bbar_remapnumbers(f, 0, 0.5) + i, xi, cmap(i,:), ...
                'FaceAlpha', p.Results.alpha, 'EdgeColor', 'none');
        end
    end
    if p.Results.showMean
        tempmean = cellfun(@(x) mean(x, 'omitnan'), rawData);
        plot((1:nA) - p.Results.BoxWidth/2, tempmean, 'k', 'LineWidth', 1.5);
    end
    if p.Results.showMedian
        tempmedian = cellfun(@(x) median(x, 'omitnan'), rawData);
        plot((1:nA) - p.Results.BoxWidth/2, tempmedian, 'k:', 'LineWidth', 1.5);
    end

    ax1 = gca;
    if isempty(p.Results.StatsYPos)
        distY = ax1.YTick(end) + (ax1.YTick(end) - ax1.YTick(end-1)) * (1:nPairs)';
    else
        distY = zeros(nPairs, 1);
        for i = 1:nPairs
            if i == 1; distY(i) = p.Results.StatsYPos(1);
            else;       distY(i) = distY(i-1) + p.Results.StatsYPos(2);
            end
        end
    end

    for i = 1:nPairs
        pAdj = PValue(i);
        if     pAdj <= 0.0001; pStr = '****';
        elseif pAdj <= 0.001;  pStr = '***';
        elseif pAdj <= 0.01;   pStr = '**';
        elseif pAdj <= 0.05;   pStr = '*';
        elseif pAdj <= 0.1;    pStr = '·';
        else;                  pStr = 'ns';
        end
        if p.Results.printP
            if round(pAdj, 4) == 0
                printedPval = ['p < 0.0001 ' pStr];
            else
                printedPval = ['p = ' num2str(round(pAdj, 4), '%.4f') ' ' pStr];
            end
        else
            printedPval = pStr;
        end
        x1c = compMatrix(i, 1); x2c = compMatrix(i, 2);
        plot([x1c x2c], repmat(distY(i), 1, 2), 'k');
        text((x1c + x2c)/2, distY(i), printedPval, ...
            'HorizontalAlignment', 'Center', 'VerticalAlignment', 'top', ...
            'FontName', p.Results.fontName, 'FontSize', p.Results.fontSize);
    end

    xticks(1:nA);
    if ~isempty(p.Results.xtickslabels) && length(p.Results.xtickslabels) == nA
        xticklabels(p.Results.xtickslabels);
    else
        if ~isempty(p.Results.xtickslabels)
            warning('Number of xtickslabels (%d) does not match number of levels (%d). Using model level names.', ...
                length(p.Results.xtickslabels), nA);
        end
        xticklabels(levelsA);
    end
    xlim([0.5, nA + 0.75]);
    ylabel(p.Results.ylabel);
    xlabel(p.Results.xlabel);
    title(p.Results.title);
    if strcmp(p.Results.grid, 'on'); grid on; else; grid off; end
    set(gcf, 'Color', p.Results.bgColor);
    set(gca, 'FontName', p.Results.fontName, 'FontSize', p.Results.fontSize);
    return
end

% CASE 2: Categorical x Continuous — Johnson-Neyman Technique
% =========================================================

contData = model.Variables.(predB);
contRange = linspace(min(contData), max(contData), p.Results.nPoints)';
sigAlpha  = p.Results.sigAlpha;

% Reference level is levelsA{1} implicitly
% We compute the simple slope of predB for each level of predA,
% and the difference between slopes (i.e. the interaction effect)
% as a function of predB, then find JN transition points.

% Find relevant coefficient indices
idxMainA   = find(contains(coeffNames, predA) & ~contains(coeffNames, ':'));
idxMainB   = find(contains(coeffNames, predB) & ~contains(coeffNames, ':'));
idxIntAB   = find(contains(coeffNames, predA) &  contains(coeffNames, predB) & contains(coeffNames, ':'));

nPairs  = nchoosek(nA, 2);
results = struct();

% Open figure
if doPlot
    if p.Results.nexttile == 0
        figure('Position', p.Results.fPosition);
        if p.Results.doPlotJN
            tiledlayout(1, 2, 'TileSpacing', 'tight', 'Padding', 'tight');
        else
            tiledlayout(1, 1, 'TileSpacing', 'tight', 'Padding', 'tight');
        end
    end
end

% ----------------------------------------------------------
% Panel 1: Regression lines per level of predA
% ----------------------------------------------------------
if doPlot
    nexttile; hold on;
    responseData = model.Variables.(model.ResponseName);

    for i = 1:nA
        % Scatter raw data
        mask = strcmp(string(model.Variables.(predA)), levelsA{i});
        scatter(contData(mask), responseData(mask), p.Results.ScatterSize, ...
            cmap(i,:), 'filled', 'MarkerFaceAlpha', p.Results.alpha);

        % Build prediction table at mean of other variables
        predTable = buildPredTable(model, predA, levelsA{i}, predB, contRange);
        [yPred, yCI] = predict(model, predTable);

        % Plot regression line and CI
        fill([contRange; flipud(contRange)], [yCI(:,1); flipud(yCI(:,2))], ...
            cmap(i,:), 'FaceAlpha', p.Results.alpha * 0.5, 'EdgeColor', 'none');
        plot(contRange, yPred, 'Color', cmap(i,:), 'LineWidth', 2);
    end

    ylabel(p.Results.ylabel);
    xlabel(predB);
    title([predA ' x ' predB]);
    if strcmp(p.Results.grid, 'on'); grid on; else; grid off; end
    set(gca, 'FontName', p.Results.fontName, 'FontSize', p.Results.fontSize);
end

% ----------------------------------------------------------
% Panel 2: Johnson-Neyman — where is the difference significant?
% ----------------------------------------------------------

% For each pair of levels of predA, compute the difference in predicted
% values across the range of predB, with SE, t, and p at each point.

pairCount = 0;
allJN = struct();

for i = 1:nA
    for j = i+1:nA
        pairCount = pairCount + 1;

        tVals  = zeros(p.Results.nPoints, 1);
        pVals  = zeros(p.Results.nPoints, 1);
        diffVals = zeros(p.Results.nPoints, 1);
        seVals = zeros(p.Results.nPoints, 1);

        for pt = 1:p.Results.nPoints
            xval = contRange(pt);

            % Contrast vector: difference between level i and j of predA
            % at this specific value of predB
            % diff = (b_Ai - b_Aj) + (b_Ai:B - b_Aj:B) * xval
            C = zeros(1, model.NumCoefficients);

            % Main effect of predA
            if i > 1; C(idxMainA(i-1)) =  1; end
            if j > 1; C(idxMainA(j-1)) = -1; end

            % Interaction term scaled by xval
            intIdx_i = findIntIdx(coeffNames, predA, levelsA{i}, predB);
            intIdx_j = findIntIdx(coeffNames, predA, levelsA{j}, predB);
            if ~isempty(intIdx_i); C(intIdx_i) = C(intIdx_i) + xval; end
            if ~isempty(intIdx_j); C(intIdx_j) = C(intIdx_j) - xval; end

            diffVals(pt) = C * fe;
            seVals(pt)   = sqrt(max(C * model.CoefficientCovariance * C', 0));
            tVals(pt)    = diffVals(pt) / seVals(pt);

            % Degrees of freedom from coefTest at this point
            [pVals(pt), ~, ~, df2] = coefTest(model, C);
        end

        % Cohen's d across the range
        dVals = diffVals / sigma;

        % Johnson-Neyman: find regions where p < sigAlpha
        isSig    = pVals < sigAlpha;
        jnBounds = findJNBounds(contRange, isSig);

        % Store results
        fname = matlab.lang.makeValidName([levelsA{i} '_vs_' levelsA{j}]);
        allJN.(fname).contRange  = contRange;
        allJN.(fname).diffVals   = diffVals;
        allJN.(fname).seVals     = seVals;
        allJN.(fname).tVals      = tVals;
        allJN.(fname).pVals      = pVals;
        allJN.(fname).dVals      = dVals;
        allJN.(fname).isSig      = isSig;
        allJN.(fname).jnBounds   = jnBounds;
        allJN.(fname).df         = df2;

        fprintf('\n=== Johnson-Neyman: %s vs %s ===\n', levelsA{i}, levelsA{j});
        if isempty(jnBounds)
            fprintf('No significant region found across the range of %s.\n', predB);
        else
            for b = 1:size(jnBounds, 1)
                fprintf('Significant for %s in [%.3f, %.3f]\n', predB, jnBounds(b,1), jnBounds(b,2));
            end
        end
        fprintf('Cohen''s d range: [%.3f, %.3f]\n', min(dVals), max(dVals));
    end
end

results = allJN;

% ----------------------------------------------------------
% JN Plot: difference with SE band, significance shading
% ----------------------------------------------------------
if p.Results.doPlotJN
    nexttile; hold on;

    pairCount = 0;
    fnames = fieldnames(allJN);

    for f = 1:length(fnames)
        pairCount = pairCount + 1;
        jn = allJN.(fnames{f});

        col = cmap(min(pairCount, size(cmap,1)), :);

        % Shade significant regions
        sigRegions = jn.isSig;
        shadedX = jn.contRange;
        shadedY = jn.diffVals;
        seY     = jn.seVals;

        % Full SE band (light)
        fill([shadedX; flipud(shadedX)], ...
             [shadedY + seY; flipud(shadedY - seY)], ...
             col, 'FaceAlpha', 0.15, 'EdgeColor', 'none');

        % Significant SE band (darker)
        sigShadedX = shadedX; sigShadedX(~sigRegions) = NaN;
        fill([sigShadedX; flipud(sigShadedX)], ...
             [shadedY + seY; flipud(shadedY - seY)], ...
             col, 'FaceAlpha', 0.35, 'EdgeColor', 'none');

        % Difference line
        plot(shadedX, shadedY, 'Color', col, 'LineWidth', 2, ...
            'DisplayName', strrep(fnames{f}, '_', ' '));

        % Mark JN boundaries
        if ~isempty(jn.jnBounds)
            for b = 1:size(jn.jnBounds, 1)
                xline(jn.jnBounds(b,1), '--', 'Color', col, 'LineWidth', 1, ...
                    'Label', sprintf('%.1f', jn.jnBounds(b,1)));
                xline(jn.jnBounds(b,2), '--', 'Color', col, 'LineWidth', 1, ...
                    'Label', sprintf('%.1f', jn.jnBounds(b,2)));
            end
        end
    end

    % Zero reference line
    yline(0, 'k--', 'LineWidth', 1);

    ylabel('Difference in predicted outcome');
    xlabel(predB);
    title('Johnson-Neyman regions');
    %legend('Location', 'best');
    if strcmp(p.Results.grid, 'on'); grid on; else; grid off; end
    set(gca, 'FontName', p.Results.fontName, 'FontSize', p.Results.fontSize);
end

if doPlot
    % sgtitle([predA ' x ' predB], 'FontName', p.Results.fontName);
    set(gcf, 'Color', p.Results.bgColor);

    pairCount = 0;
    fnames = fieldnames(allJN);

    for f = 1:length(fnames)
        pairCount = pairCount + 1;
        jn = allJN.(fnames{f});

        % Mark JN boundaries
        if ~isempty(jn.jnBounds)
            for b = 1:size(jn.jnBounds, 1)
                xline(jn.jnBounds(b,1), ':', 'Color', cmap(1,:), 'LineWidth', 1, ...
                    'Label', sprintf('%.1f', jn.jnBounds(b,1)));
                xline(jn.jnBounds(b,2), ':', 'Color', cmap(1,:), 'LineWidth', 1, ...
                    'Label', sprintf('%.1f', jn.jnBounds(b,2)));
            end
        end
    end

    if ~isempty(p.Results.xtickslabels) && length(p.Results.xtickslabels) == length(levelsA)
        tmp = repmat("",length(p.Results.xtickslabels)*2+1,1); counter = 1;
        for i2 = 1:3:length(p.Results.xtickslabels)*2+1
            tmp(i2) = string(p.Results.xtickslabels{counter});
            counter = counter + 1;
        end
        legend(tmp, 'Location', 'best');
    else
        legend(levelsA, 'Location', 'best');
    end

end

end

% =========================================================
% Helper: build prediction table for one level of predA
% =========================================================
function predTable = buildPredTable(model, predA, levelA, predB, contRange)
    nPoints  = length(contRange);
    varNames = model.VariableNames;
    predTable = table();

    for v = 1:length(varNames)
        vname = varNames{v};
        if strcmp(vname, model.ResponseName); continue; end

        colData = model.Variables.(vname);

        if strcmp(vname, predA)
            if iscategorical(colData)
                predTable.(vname) = repmat(categorical({levelA}), nPoints, 1);
            else
                predTable.(vname) = repmat(string(levelA), nPoints, 1);
            end
        elseif strcmp(vname, predB)
            predTable.(vname) = contRange;
        elseif isnumeric(colData)
            predTable.(vname) = repmat(mean(colData, 'omitnan'), nPoints, 1);
        elseif iscategorical(colData)
            predTable.(vname) = repmat(mode(colData), nPoints, 1);
        else
            predTable.(vname) = repmat(colData(1), nPoints, 1);
        end
    end
end

% =========================================================
% Helper: find interaction coefficient index for predA level and predB
% =========================================================
function idx = findIntIdx(coeffNames, predA, levelA, predB)
    matchLevel = contains(coeffNames, levelA);
    matchB     = contains(coeffNames, predB);
    matchColon = contains(coeffNames, ':');
    idx = find(matchLevel & matchB & matchColon);
end

% =========================================================
% Helper: find Johnson-Neyman boundary points
% =========================================================
function bounds = findJNBounds(x, isSig)
    bounds = [];
    transitions = diff([0; isSig; 0]);
    starts = find(transitions ==  1);
    ends   = find(transitions == -1) - 1;
    for k = 1:length(starts)
        bounds(k, :) = [x(starts(k)), x(ends(k))]; %#ok<AGROW>
    end
end