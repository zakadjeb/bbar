function [allStats, windowMeans, rcDataByWindow, tl] = bbar_timeseries(data, param, method, varargin)
% BBAR_TIMESERIES
%
% INPUT
%   data   : 1 x N cell
%            N = number of conditions
%            each cell contains a matrix of size:
%               participants x time
%
%   param  : parameter passed to bbar_subdivide()
%            - window size       if method = 'windowSize'
%            - number of outputs if method = 'nSamples'
%
%   method : 'windowSize' or 'nSamples'
%
% WRAPPER-SPECIFIC NAME-VALUE PAIRS
%   'Time'             : time vector for x-axis (default = 1:T)
%   'ConditionLabels'  : 1 x N cell labels
%   'FigureTitle'      : overall tiledlayout title
%   'TopYLabel'        : y-label of top plot
%   'TopXLabel'        : x-label of top plot
%   'TopTitle'         : title of top plot
%   'TopLineWidth'     : line width of mean traces
%   'PatchAlpha'       : transparency of SD patch
%   'MarkWindows'      : draw vertical boundaries (default = true)
%   'TopLegend'        : show legend on top plot (default = true)
%
% All remaining varargin are passed directly to bbar_raincloud(), so you
% can still use:
%   'showMean', 'Scatter', 'BoxPlot', 'Normal', 'paired', 'stats',
%   'comparison', 'cmap', 'ylabel', etc.
%
% OUTPUT
%   windowMeans    : participants x conditions x windows
%   rcDataByWindow : 1 x windows cell, each cell is 1 x conditions cell
%   tl             : tiledlayout handle
%
% EXAMPLE
%   [wm, rc, tl] = bbar_timeseries( ...
%       data, 5, 'nSamples', ...
%       'ConditionLabels', {'A','B','C'}, ...
%       'FigureTitle', 'Windowed Means', ...
%       'showMean', 1, 'Scatter', 1, 'BoxPlot', 1, 'Normal', 1, ...
%       'paired', 1, 'ylabel', 'Window mean');

    %% ----------------------------
    % Validate main inputs
    % ----------------------------
    if ~iscell(data) || isempty(data) || size(data,1) ~= 1
        error('Input "data" must be a non-empty 1 x N cell array.');
    end

    if ~isscalar(param) || ~isnumeric(param) || param <= 0 || round(param) ~= param
        error('Input "param" must be a positive integer scalar.');
    end

    if ~(ischar(method) || isstring(method))
        error('Input "method" must be ''windowSize'' or ''nSamples''.');
    end
    method = char(method);

    nConditions = numel(data);

    allStats = cell(1,param);

    %% ----------------------------
    % Wrapper defaults
    % ----------------------------
    wrapper.Time            = [];
    wrapper.ConditionLabels = arrayfun(@(i) sprintf('Condition %d', i), 1:nConditions, 'UniformOutput', false);
    wrapper.FigureTitle     = '';
    wrapper.TopYLabel       = 'Amplitude';
    wrapper.TopXLabel       = 'Samples';
    wrapper.TopTitle        = 'Time series';
    wrapper.TopLineWidth    = 2;
    wrapper.PatchAlpha      = 0.18;
    wrapper.MarkWindows     = true;
    wrapper.TopLegend       = true;

    % Everything else goes to bbar_raincloud
    rainArgs = varargin;

    recognized = {'Time','ConditionLabels','FigureTitle','TopYLabel','TopXLabel', ...
                  'TopTitle','TopLineWidth','PatchAlpha','MarkWindows','TopLegend'};

    i = 1;
    while i <= numel(rainArgs)
        if ischar(rainArgs{i}) || isstring(rainArgs{i})
            key = char(rainArgs{i});
            idx = find(strcmpi(key, recognized), 1);

            if ~isempty(idx)
                if i == numel(rainArgs)
                    error('Missing value for wrapper argument "%s".', key);
                end
                val = rainArgs{i+1};

                switch lower(key)
                    case 'time'
                        wrapper.Time = val;
                    case 'conditionlabels'
                        wrapper.ConditionLabels = val;
                    case 'figuretitle'
                        wrapper.FigureTitle = val;
                    case 'topylabel'
                        wrapper.TopYLabel = val;
                    case 'topxlabel'
                        wrapper.TopXLabel = val;
                    case 'toptitle'
                        wrapper.TopTitle = val;
                    case 'toplinewidth'
                        wrapper.TopLineWidth = val;
                    case 'patchalpha'
                        wrapper.PatchAlpha = val;
                    case 'markwindows'
                        wrapper.MarkWindows = val;
                    case 'toplegend'
                        wrapper.TopLegend = val;
                end

                rainArgs(i:i+1) = [];
                continue;
            end
        end
        i = i + 1;
    end

    if ~iscell(wrapper.ConditionLabels) || numel(wrapper.ConditionLabels) ~= nConditions
        error('ConditionLabels must be a 1 x N cell array matching the number of conditions.');
    end

    %% ----------------------------
    % Validate condition matrices
    % ----------------------------
    nParticipants = [];
    tsLength = [];

    for c = 1:nConditions
        x = data{c};

        if ~isnumeric(x) || isempty(x) || ndims(x) ~= 2
            error('Each data{c} must contain a non-empty numeric matrix of size participants x time.');
        end

        [p, t] = size(x);

        if isempty(nParticipants)
            nParticipants = p;
            tsLength = t;
        else
            if p ~= nParticipants || t ~= tsLength
                error('All condition matrices must have the same size: participants x time.');
            end
        end
    end

    if isempty(wrapper.Time)
        wrapper.Time = 1:tsLength;
    else
        if numel(wrapper.Time) ~= tsLength
            error('The Time vector must match the number of time samples.');
        end
        wrapper.Time = wrapper.Time(:)';
    end

    %% ----------------------------
    % Use bbar_subdivide() on first participant to determine W
    % ----------------------------
    testSub = bbar_subdivide(data{1}(1,:), param, method);
    nWindows = numel(testSub);

    % windowMeans(participant, condition, window)
    windowMeans = nan(nParticipants, nConditions, nWindows);

    for c = 1:nConditions
        for p = 1:nParticipants
            tmp = bbar_subdivide(data{c}(p,:), param, method);

            if numel(tmp) ~= nWindows
                error(['bbar_subdivide() returned different numbers of windows across inputs. ' ...
                       'Make sure all time series have identical length.']);
            end

            windowMeans(p,c,:) = tmp(:);
        end
    end

    %% ----------------------------
    % Reconstruct boundaries for top-plot window markers
    % (same logic as bbar_subdivide)
    % ----------------------------
    [starts, ends] = i_getWindowEdges(tsLength, param, method);
    if numel(starts) ~= nWindows
        error('Internal mismatch: reconstructed window count does not match bbar_subdivide().');
    end

    %% ----------------------------
    % Prepare raincloud input by window
    % ----------------------------
    rcDataByWindow = cell(1, nWindows);
    for w = 1:nWindows
        rcData = cell(1, nConditions);
        for c = 1:nConditions
            rcData{c} = squeeze(windowMeans(:,c,w));
        end
        rcDataByWindow{w} = rcData;
    end

    %% ----------------------------
    % Determine colormap if user passed one to bbar_raincloud
    % ----------------------------
    cmapUser = [];
    cmapIdx = find_namedarg(rainArgs, 'cmap');
    if ~isempty(cmapIdx)
        cmapUser = rainArgs{cmapIdx+1};
        if size(cmapUser,1) < nConditions
            error('Provided cmap must have at least as many rows as conditions.');
        end
    else
        % fallback similar spirit to your raincloud function
        cmapUser = lines(nConditions);
    end

    %% ----------------------------
    % Figure + tiledlayout
    % ----------------------------
    figure;
    tl = tiledlayout(2, nWindows, 'TileSpacing', 'compact', 'Padding', 'compact');

    %% ----------------------------
    % Top plot spanning all columns
    % ----------------------------
    axTop = nexttile(tl, 1, [1 nWindows]);
    hold(axTop, 'on');

    for c = 1:nConditions
        mu = mean(data{c}, 1, 'omitnan');
        sd = std(data{c}, 0, 1, 'omitnan');

        x = wrapper.Time(:)';
        yUpper = mu + sd;
        yLower = mu - sd;

        % Transparent SD patch
        patch(axTop, ...
            [x, fliplr(x)], ...
            [yUpper, fliplr(yLower)], ...
            cmapUser(c,:), ...
            'FaceAlpha', wrapper.PatchAlpha, ...
            'EdgeColor', 'none', ...
            'HandleVisibility', 'off');

        % Mean line
        plot(axTop, x, mu, ...
            'Color', cmapUser(c,:), ...
            'LineWidth', wrapper.TopLineWidth, ...
            'DisplayName', wrapper.ConditionLabels{c});
    end

    % Window markers
    if wrapper.MarkWindows && nWindows > 1
        for w = 1:(nWindows-1)
            xline(axTop, wrapper.Time(ends(w)), ':', 'HandleVisibility', 'off');
        end
    end

    xlabel(axTop, wrapper.TopXLabel);
    ylabel(axTop, wrapper.TopYLabel);
    title(axTop, wrapper.TopTitle);
    grid(axTop, 'off');

    if wrapper.TopLegend
        legend(axTop, 'Location', 'best');
    end

    hold(axTop, 'off');

    if ~isempty(wrapper.FigureTitle)
        title(tl, wrapper.FigureTitle, 'FontWeight', 'bold');
    end

    %% ----------------------------
    % Bottom rainclouds
    % ----------------------------
    for w = 1:nWindows
        rcArgsThis = rainArgs;

        % Add defaults only if user did not already specify them
        if ~has_namedarg(rcArgsThis, 'nexttile')
            rcArgsThis = [rcArgsThis, {'nexttile', nWindows + w}];
        end
        if ~has_namedarg(rcArgsThis, 'title')
            rcArgsThis = [rcArgsThis, {'title', sprintf('Window %d', w)}];
        end
        if ~has_namedarg(rcArgsThis, 'xtickslabels')
            rcArgsThis = [rcArgsThis, {'xtickslabels', [wrapper.ConditionLabels, {''}]}];
        end
        if ~has_namedarg(rcArgsThis, 'cmap')
            rcArgsThis = [rcArgsThis, {'cmap', cmapUser(1:nConditions,:)}];
        end

        [~, stats] = bbar_raincloud(rcDataByWindow{w}, rcArgsThis{:});
        allStats{w} = stats;
    end
end


%% =========================================================
% Helper: reconstruct same edges as bbar_subdivide()
%% =========================================================
function [starts, ends] = i_getWindowEdges(L, param, method)

    switch lower(method)
        case 'windowsize'
            winSize = param;

            if winSize >= L
                starts = 1;
                ends   = L;
                return;
            end

            starts = 1:winSize:L;
            ends   = [starts(2:end)-1, L];

            if numel(starts) > 1
                lastLen = ends(end) - starts(end) + 1;
                if lastLen < winSize
                    ends(end-1) = ends(end);
                    starts(end) = [];
                    ends(end) = [];
                end
            end

        case 'nsamples'
            nSamples = param;

            if nSamples >= L
                starts = 1:L;
                ends   = 1:L;
                return;
            end

            edges = round(linspace(1, L+1, nSamples+1));
            starts = edges(1:end-1);
            ends   = edges(2:end)-1;

            valid = ends >= starts;
            starts = starts(valid);
            ends   = ends(valid);

            if numel(starts) > 1
                winLens = ends - starts + 1;
                targetLen = round(median(winLens));

                lastLen = winLens(end);
                if lastLen < targetLen
                    ends(end-1) = ends(end);
                    starts(end) = [];
                    ends(end) = [];
                end
            end

        otherwise
            error('method must be ''windowSize'' or ''nSamples''.');
    end
end


%% =========================================================
% Helper: does name-value arg exist?
%% =========================================================
function tf = has_namedarg(args, name)
    tf = false;
    for k = 1:numel(args)-1
        if (ischar(args{k}) || isstring(args{k})) && strcmpi(char(args{k}), name)
            tf = true;
            return;
        end
    end
end


%% =========================================================
% Helper: find name-value arg index
%% =========================================================
function idx = find_namedarg(args, name)
    idx = [];
    for k = 1:numel(args)-1
        if (ischar(args{k}) || isstring(args{k})) && strcmpi(char(args{k}), name)
            idx = k;
            return;
        end
    end
end