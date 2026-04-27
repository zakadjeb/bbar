function [ax, imh, xgrid, ygrid, F] = bbar_heatmap(X, Y, varargin)
%KDE2D_HEATMAP  Plot a 2-D kernel density heatmap from X,Y samples.
%
% [ax, imh, xgrid, ygrid, F] = kde2d_heatmap(X, Y, Name,Value,...)
%
% Required:
%   X, Y        - equal-length vectors of samples
%
% Name-Value parameters (parser-validated):
%   'Resolution'     (scalar)     grid size per axis (default 200)
%   'XLim'           (2-el vec)   [xmin xmax] for grid (default data range)
%   'YLim'           (2-el vec)   [ymin ymax] for grid (default data range)
%   'Scale'          (scalar)     Blobscale (e.g., 0.5, default = auto)
%   'Kernel'         (char)       'normal'|'epanechnikov'|'box'|'triangle'
%   'Support'        (char|2x2)   'unbounded'|'positive'|'bounded'
%   'Weights'        (vector)     observation weights for ksdensity
%   'Parent'         (axes)       target axes; if empty, chooses below
%   'Tile'           (scalar)     use nexttile(Tile) if Parent empty
%   'ShowColorbar'   (logical)    show colorbar (default true)
%   'CLim'           (1x2)        color limits for imagesc
%   'Colormap'       (colormap)   e.g., parula(256) (default: parula)
%   'Title'          (char)       plot title (default '')
%   'XLabel'         (char)       x-label (default 'X')
%   'YLabel'         (char)       y-label (default 'Y')
%
% Returns:
%   ax     - axes handle used
%   imh    - image handle from imagesc
%   xgrid, ygrid - 1D grid vectors
%   F      - density matrix (rows = y, cols = x)

% -------- input checks & parser --------
validateattributes(X,{'numeric'},{'vector','nonempty'},mfilename,'X',1);
validateattributes(Y,{'numeric'},{'vector','nonempty'},mfilename,'Y',2);
assert(numel(X)==numel(Y),'X and Y must be the same length.');
X = X(:); Y = Y(:);

p = inputParser; p.FunctionName = mfilename;

addParameter(p,'Plot',true,@islogical);
addParameter(p,'Resolution',200,@(v)isnumeric(v)&&isscalar(v)&&v>=10);
addParameter(p,'XLim',[],@(v)isnumeric(v)&&isempty(v)|| (isvector(v)&&numel(v)==2));
addParameter(p,'YLim',[],@(v)isnumeric(v)&&isempty(v)|| (isvector(v)&&numel(v)==2));
addParameter(p,'Scale',[],@(v)isempty(v) || (isscalar(v)&&isnumeric(v)));
addParameter(p,'Kernel','normal',@(s)ischar(s)||isstring(s));
addParameter(p,'Support','unbounded',@(v)ischar(v)||isstring(v)|| ...
    (isnumeric(v)&&isequal(size(v),[2 2])));
addParameter(p,'Weights',[],@(v)isnumeric(v)&& (isempty(v)||isvector(v)&&numel(v)==numel(X)));

addParameter(p,'Parent',[],@(h) isempty(h) || isgraphics(h,'axes'));
addParameter(p,'Tile',[],@(v) isempty(v) || (isscalar(v)&&isnumeric(v)));

addParameter(p,'ShowColorbar',true,@islogical);
addParameter(p,'CLim',[],@(v)isnumeric(v)&& (isempty(v)|| (isvector(v)&&numel(v)==2)));
addParameter(p,'Colormap',parula(256),@(v)isnumeric(v)&&size(v,2)==3);
addParameter(p,'Title','',@(s)ischar(s)||isstring(s));
addParameter(p,'XLabel','X',@(s)ischar(s)||isstring(s));
addParameter(p,'YLabel','Y',@(s)ischar(s)||isstring(s));

parse(p,varargin{:});
S = p.Results;

% -------- grid --------
if isempty(S.XLim), S.XLim = [min(X) max(X)]; end
if isempty(S.YLim), S.YLim = [min(Y) max(Y)]; end
nx = max(10, round(S.Resolution));
ny = nx;

xgrid = linspace(S.XLim(1), S.XLim(2), nx);
ygrid = linspace(S.YLim(1), S.YLim(2), ny);
[Xg,Yg] = meshgrid(xgrid, ygrid);

% -------- ksdensity options assembly --------
pts = [X Y];
q   = [Xg(:) Yg(:)];
kdOpts = {'Kernel', char(S.Kernel), 'Support', S.Support};

if ~isempty(S.Scale)
    kdOpts = [kdOpts, {'Bandwidth', [S.Scale*std(X), S.Scale*std(Y)]}]; 
end
if ~isempty(S.Weights)
    kdOpts = [kdOpts, {'Weights', S.Weights(:)}]; 
end

% compute density on grid
f = ksdensity(pts, q, kdOpts{:});   % f is column vector
F = reshape(f, ny, nx);             % rows = y, cols = x

% -------- Hadnling plot --------
if S.Plot
    % -------- target axes selection --------
    if ~isempty(S.Parent)
        ax = S.Parent;
    elseif ~isempty(S.Tile)
        ax = nexttile(S.Tile);
    else
        figure('Color','w');
        ax = axes('Box','on');
    end
    
    % -------- plot --------
    imh = imagesc(ax, xgrid, ygrid, F);
    axis(ax, 'xy');  % keep y increasing upward
    xlabel(ax, S.XLabel); ylabel(ax, S.YLabel);
    if ~isempty(S.Title), title(ax, S.Title); end
    colormap(ax, S.Colormap);
    if ~isempty(S.CLim), set(ax,'CLim',S.CLim); end
    if S.ShowColorbar, colorbar(ax); end
else
    ax = NaN;
    imh = NaN;
end

end