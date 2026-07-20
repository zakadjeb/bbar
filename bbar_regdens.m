function bbar_regdens(data, varargin)
% BBAR_REGDENS  Density heatmap (or contour) with scatter and optional
%               regression overlay — all in one axes.
%
%   bbar_regdens(data)
%   bbar_regdens(data, Name, Value, ...)
%
%   INPUT
%     data   : cell array, where each cell is an n x 2 numeric matrix.
%              Column 1 = x values, Column 2 = y values.
%
%              Example:
%                data{1} = [2,2; 2,3; 1,3];
%                data{2} = [1,5; 2,1; 3,3];
%
%   OPTIONAL NAME-VALUE PAIRS
%     'yLim'          : [min max]  y-axis limits            (default: auto)
%     'xLim'          : [min max]  x-axis limits            (default: auto)
%     'resolution'    : scalar     KDE grid resolution      (default: 300)
%     'scale'         : scalar     kernel bandwidth scale   (default: 0.3)
%     'xlabel'        : string     x-axis label             (default: 'x')
%     'ylabel'        : string     y-axis label             (default: 'y')
%     'title'         : string     plot title               (default: '')
%     'ShowColorbar'  : logical    show colorbar            (default: false)
%     'colormap'      : string     slanCM colormap name     (default: 'BuPu')
%     'MarkerSize'    : scalar     scatter marker size      (default: 6)
%     'MarkerAlpha'   : scalar     scatter marker alpha     (default: 0.5)
%     'ShowRegression': logical    OLS line + 95 % CI       (default: false)
%     'nPred'         : scalar     x points for regression  (default: 200)
%     'density'       : 'heatmap' | 'contour'              (default: 'heatmap')
%     'nContours'     : scalar     contour levels           (default: 3)
%     'tile'          : scalar | [] nexttile index, [] = new figure (default: [])
%     'DensityScatter': logical    colour each marker by its local KDE density,
%                                  using the same colormap as the background.
%                                  Low-density points → dark end of cmap;
%                                  high-density points → bright end. (default: false)
%     'ShowDensity'   : logical    show the density layer (heatmap or contour).
%                                  Set to false for scatter/regression only, while
%                                  still using the KDE for DensityScatter. (default: true)
%
%   LAYER ORDER (bottom to top)
%     1. Density heatmap or contour   [from bbar_heatmap ax]
%     2. 95 % CI patch per cell
%     3. Regression line per cell
%     4. Scatter markers per cell

% -------------------------------------------------------------------------
% 1. Parse inputs
% -------------------------------------------------------------------------
p = inputParser;
p.KeepUnmatched = true;

validData = @(c) iscell(c) && all(cellfun(@(m) ...
    isnumeric(m) && ismatrix(m) && size(m,2)==2, c));
addRequired(p, 'data', validData);

addParameter(p, 'yLim',          [],        @(x) isempty(x)||(isnumeric(x)&&numel(x)==2));
addParameter(p, 'xLim',          [],        @(x) isempty(x)||(isnumeric(x)&&numel(x)==2));
addParameter(p, 'resolution',    300,       @isnumeric);
addParameter(p, 'scale',         0.3,       @isnumeric);
addParameter(p, 'ShowColorbar',  false,     @islogical);
addParameter(p, 'xlabel',        'x',       @ischar);
addParameter(p, 'ylabel',        'y',       @ischar);
addParameter(p, 'title',         '',        @ischar);
addParameter(p, 'colormap',      'BuPu',    @ischar);
addParameter(p, 'MarkerSize',    3,         @isnumeric);
addParameter(p, 'MarkerAlpha',   .5,        @isnumeric);
addParameter(p, 'ShowRegression',false,     @islogical);
addParameter(p, 'nPred',         200,       @isnumeric);
addParameter(p, 'density',       'heatmap', @(s) ismember(s,{'heatmap','contour'}));
addParameter(p, 'nContours',     3,         @isnumeric);
addParameter(p, 'tile',          [],        @(x) isempty(x)||isnumeric(x));
addParameter(p, 'DensityScatter',false,     @islogical);
addParameter(p, 'ShowDensity',   true,      @islogical);

parse(p, data, varargin{:});
opts   = p.Results;
nCells = numel(data);

% -------------------------------------------------------------------------
% 2. Pool all data + auto limits
% -------------------------------------------------------------------------
xAll = cell2mat(cellfun(@(m) m(:,1), data(:), 'UniformOutput', false));
yAll = cell2mat(cellfun(@(m) m(:,2), data(:), 'UniformOutput', false));

padFrac = 0.05;
if isempty(opts.xLim)
    r = max(xAll) - min(xAll);
    opts.xLim = [min(xAll) - padFrac*r,  max(xAll) + padFrac*r];
end
if isempty(opts.yLim)
    r = max(yAll) - min(yAll);
    opts.yLim = [min(yAll) - padFrac*r,  max(yAll) + padFrac*r];
end

% -------------------------------------------------------------------------
% 3. Colormap + per-cell colours
% -------------------------------------------------------------------------
cmap    = slanCM(opts.colormap);
fgColor = abs(1 - cmap(1,:));

cellIdx  = round(linspace(round(size(cmap,1)*0.4), size(cmap,1), nCells));
cellCols = cmap(cellIdx, :);

% -------------------------------------------------------------------------
% 4. Layer 1 — density base
%    bbar_heatmap returns [ax, imh, xgrid, ygrid, F]
%    We use its ax handle for all subsequent drawing.
% -------------------------------------------------------------------------
if ~isempty(opts.tile)
    tileArg = {'tile', opts.tile};
else
    tileArg = {};
end

switch opts.density

    case 'heatmap'
        [ax, imh, xgrid, ygrid, F] = bbar_heatmap(xAll, yAll, ...
            'yLim',         opts.yLim, ...
            'xLim',         opts.xLim, ...
            'resolution',   opts.resolution, ...
            'scale',        opts.scale, ...
            'ylabel',       opts.ylabel, ...
            'xlabel',       opts.xlabel, ...
            'ShowColorbar', opts.ShowColorbar && opts.ShowDensity, ...
            tileArg{:});

        % Hide the image layer if ShowDensity is off (KDE grid still available)
        if ~opts.ShowDensity
            imh.Visible = 'off';
        end

    case 'contour'
        [ax, ~, xgrid, ygrid, F] = bbar_heatmap(xAll, yAll, ...
            'yLim',         opts.yLim, ...
            'xLim',         opts.xLim, ...
            'resolution',   opts.resolution, ...
            'scale',        opts.scale, ...
            'ylabel',       opts.ylabel, ...
            'xlabel',       opts.xlabel, ...
            'ShowColorbar', false, ...
            tileArg{:});

        % Replace image with contours drawn on the same ax
        cla(ax);
        if opts.ShowDensity
            contourf(ax, xgrid, ygrid, F, opts.nContours, 'LineColor','none');
            contour( ax, xgrid, ygrid, F, opts.nContours, ...
                     'LineColor', abs(1-cmap(1,:)), 'LineWidth', 0.8);
            if opts.ShowColorbar
                colorbar(ax);
            end
        end
        xlim(ax, opts.xLim);
        ylim(ax, opts.yLim);
end

% Apply colormap and figure background colour
f = ancestor(ax, 'figure');
colormap(ax, cmap);
f.Color = [1 1 1];

% -------------------------------------------------------------------------
% 5. Layers 2-4 — per cell: CI patch → regression line → scatter
%    All drawn onto the ax returned by bbar_heatmap.
% -------------------------------------------------------------------------
hold(ax, 'on');
densPerCell = cell(nCells, 1);   % pre-allocate for DensityScatter

for k = 1:nCells
    x = data{k}(:,1);
    y = data{k}(:,2);
    c = cellCols(k,:);

    if opts.ShowRegression
        [xPred, yPred, yCI] = ols_with_ci(x, y, opts.xLim, opts.nPred);

        % Layer 2 — CI patch
        patch(ax, [xPred; flipud(xPred)], ...
                  [yCI(:,1); flipud(yCI(:,2))], ...
                  c, 'FaceAlpha', 0.25, 'EdgeColor', 'none');

        % Layer 3 — regression line
        plot(ax, xPred, yPred, '-', 'LineWidth', 2, 'Color', c);
    end

    if opts.DensityScatter
        % Sample KDE at each point's location; draw after loop for global normalisation
        densPerCell{k} = interp2(xgrid, ygrid, F, x, y, 'linear', 0);
    else
        % Layer 4 — standard scatter
        scatter(ax, x, y, opts.MarkerSize, cmap(end,:), 'filled', ...
            'MarkerFaceAlpha', opts.MarkerAlpha, 'MarkerEdgeColor', 'none');
    end
end

% -------------------------------------------------------------------------
% 5b. Density-coloured scatter
%     Drawn after the loop so colour normalisation is global across all cells.
% -------------------------------------------------------------------------
if opts.DensityScatter
    allDens = cell2mat(densPerCell);
    dMin    = min(allDens);
    dMax    = max(allDens);
    nColors = size(cmap, 1);

    for k = 1:nCells
        x  = data{k}(:,1);
        y  = data{k}(:,2);
        pd = densPerCell{k};

        % Map each point's density to a colormap row
        normD  = (pd - dMin) / max(dMax - dMin, eps);    % [0, 1]
        cidx   = max(1, round(normD * (nColors-1)) + 1); % [1, nColors]
        pcolor = cmap(cidx, :);                           % n x 3 RGB

        % scatter() accepts per-point CData as n x 3
        scatter(ax, x, y, opts.MarkerSize, pcolor, 'filled', ...
            'MarkerFaceAlpha', opts.MarkerAlpha, 'MarkerEdgeColor', 'none');
    end
end

% -------------------------------------------------------------------------
% 6. Decoration
% -------------------------------------------------------------------------
yline(ax, 0, ':', 'Color', fgColor);
xlabel(ax, opts.xlabel);
ylabel(ax, opts.ylabel);
if ~isempty(opts.title)
    title(ax, opts.title);
end

end % main function


% =========================================================================
% LOCAL — OLS regression + 95 % CI (no toolboxes)
% =========================================================================
function [xPred, yPred, yCI] = ols_with_ci(x, y, xLim, nPred)
% y = b0 + b1*x fitted by OLS; 95 % CI on the mean prediction.
%
%   SE(x0) = s * sqrt( 1/n + (x0-xbar)^2 / Sxx )
%   CI     = yhat +/- t(0.975, n-2) * SE(x0)

n     = numel(x);
X     = [ones(n,1), x(:)];
b     = X \ y(:);

xPred = linspace(xLim(1), xLim(2), nPred)';
yPred = [ones(nPred,1), xPred] * b;

RSS    = sum((y(:) - X*b).^2);
s2     = RSS / (n - 2);
xbar   = mean(x);
Sxx    = sum((x(:) - xbar).^2);
sePred = sqrt(s2 * (1/n + (xPred - xbar).^2 / Sxx));

tCrit = 1.96;
if exist('tinv','file')
    tCrit = tinv(0.975, n - 2);
end

yCI = [yPred - tCrit.*sePred,  yPred + tCrit.*sePred];

end