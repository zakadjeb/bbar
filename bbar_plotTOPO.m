function [Stats, topoData, topoSubjectData, rcData, f] = bbar_plotTOPO(data, varargin)
% BBAR_PLOTTOPO
% Simple topoplot function for EEGLAB ALLEEG-style data.
%
% Input example:
%   data = {ALLEEG(IdxRewardEroticPre) ALLEEG(IdxRewardEroticPost)};
%
% Each cell is one condition.
% Each condition contains an array/list of EEG structs.
% Each EEG struct is treated as one participant.
%
% The function:
%   1. bandpass filters each participant if freqWindow is given
%   2. extracts a timeWindow
%   3. computes one value per channel
%   4. averages across participants per condition
%   5. plots topoplots using topoplot()
%   6. optionally extracts channel/ROI values for bbar_raincloud()
%   7. for an ROI use: 'chan',{'Cz','C1','C3'} 
%
% Outputs:
%   Stats            output from bbar_raincloud, if used
%   topoData         1 x nConditions cell; each cell is 1 x channels group mean
%   topoSubjectData  1 x nConditions cell; each cell is subjects x channels
%   rcData           raincloud data; one vector per condition
%
% Example:
%   [Stats, topoData, topoSub, rcData] = bbar_plotTOPO( ...
%       {ALLEEG(IdxRewardEroticPre) ALLEEG(IdxRewardEroticPost)}, ...
%       'freqWindow', [8 12], ...
%       'timeWindow', [300 600], ...
%       'measure', 'power', ...
%       'chan', 'Cz', ...
%       'raincloud', 1, ...
%       'names', {'Pre','Post'}, ...
%       'raincloudArgs', {'paired', 1}, ...
%       'topoplotArgs', {'style','both','electrodes','on'});

% -------------------- Input parser --------------------
p = inputParser;
addRequired(p, 'data', @iscell);

addParameter(p, 'freqWindow', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
addParameter(p, 'timeWindow', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));

addParameter(p, 'measure', 'power', @(x) ischar(x) || isstring(x));
addParameter(p, 'chan', [], @(x) isempty(x) || isnumeric(x) || ischar(x) || isstring(x) || iscell(x));

addParameter(p, 'raincloud', 0, @(x) isnumeric(x) || islogical(x));
addParameter(p, 'raincloudArgs', {}, @iscell);

addParameter(p, 'names', {}, @iscell);
addParameter(p, 'plot', 1, @(x) isnumeric(x) || islogical(x));
addParameter(p, 'newfig', 1, @(x) isnumeric(x) || islogical(x));
addParameter(p, 'maplimits', [], @(x) isempty(x) || isnumeric(x) || ischar(x) || isstring(x));
addParameter(p, 'topoplotArgs', {}, @iscell);
addParameter(p, 'filterOrder', 4, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'markChannels', 1, @(x) isnumeric(x) || islogical(x));
addParameter(p, 'markMarker', 'x', @(x) ischar(x) || isstring(x));
addParameter(p, 'markColor', 'w');
addParameter(p, 'markSize', 7, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'markLineWidth', 2, @(x) isnumeric(x) && isscalar(x));
addOptional(p, 'cmap', jet(size(data,2)),...
    @(x) (size(data,2) <= size(x,1)) && isnumeric(x));

parse(p, data, varargin{:});

freqWindow = p.Results.freqWindow;
timeWindow = p.Results.timeWindow;
measure = lower(char(p.Results.measure));

Stats = table();
rcData = {};

% -------------------- Basic info from first EEG --------------------
EEG0 = data{1}(1);

times = double(EEG0.times);
srate = double(EEG0.srate);
chanlocs = EEG0.chanlocs;
nChan = size(EEG0.data, 1);

if isempty(timeWindow)
    timeWindow = [times(1) times(end)];
end

timeIdx = times >= timeWindow(1) & times <= timeWindow(2);

if ~any(timeIdx)
    error('timeWindow does not overlap EEG.times.');
end

doFilter = ~isempty(freqWindow);

if doFilter
    if freqWindow(2) >= srate/2
        error('freqWindow upper limit must be below Nyquist frequency.');
    end

    [b, a] = butter(p.Results.filterOrder, freqWindow ./ (srate/2), 'bandpass');
end

% -------------------- Selected channels for raincloud and topo marking --------------------
chanIdx = [];

if ~isempty(p.Results.chan)
    chanIdx = find_channels(p.Results.chan, chanlocs);
end

% -------------------- Condition names --------------------
nCond = numel(data);

if isempty(p.Results.names)
    names = cell(1, nCond);
    for c = 1:nCond
        names{c} = sprintf('Condition %d', c);
    end
else
    names = p.Results.names;
end

% -------------------- Extract topographies --------------------
topoSubjectData = cell(1, nCond);
topoData = cell(1, nCond);

for c = 1:nCond

    EEGs = data{c};

    if iscell(EEGs)
        EEGs = [EEGs{:}];
    end

    nSub = numel(EEGs);
    topoSubjectData{c} = nan(nSub, nChan);

    for s = 1:nSub

        X = double(EEGs(s).data); % channels x time x trials OR channels x time

        if ndims(X) == 2
            X = reshape(X, size(X,1), size(X,2), 1);
        end

        % ----- Bandpass filter -----
        if doFilter
            nTrials = size(X, 3);

            % filtfilt works along rows, so make data time x observations
            Xr = reshape(permute(X, [2 1 3]), size(X,2), nChan*nTrials);
            Xr = filtfilt(b, a, Xr);

            % back to channels x time x trials
            X = permute(reshape(Xr, size(X,2), nChan, nTrials), [2 1 3]);
        end

        % ----- Extract time window -----
        Xwin = X(:, timeIdx, :);
        Xwin = reshape(Xwin, nChan, []); % channels x samples

        % ----- Reduce to one value per channel -----
        switch measure
            case 'mean'
                vals = mean(Xwin, 2, 'omitnan');

            case 'absmean'
                vals = mean(abs(Xwin), 2, 'omitnan');

            case 'power'
                vals = mean(Xwin.^2, 2, 'omitnan');

            case 'rms'
                vals = sqrt(mean(Xwin.^2, 2, 'omitnan'));

            otherwise
                error('Unknown measure. Use mean, absmean, power, or rms.');
        end

        topoSubjectData{c}(s,:) = vals(:)';
    end

    % Mean topography across participants
    topoData{c} = mean(topoSubjectData{c}, 1, 'omitnan');
end

% -------------------- Raincloud data --------------------
if p.Results.raincloud

    if isempty(chanIdx)
        error('To use raincloud=1, specify a channel or ROI using ''chan'', e.g. ''Cz'' or {''Cz'',''Pz''}.');
    end

    rcData = cell(1, nCond);

    for c = 1:nCond
        % one value per participant
        % if multiple channels are given, average across those channels
        rcData{c} = mean(topoSubjectData{c}(:, chanIdx), 2, 'omitnan');
    end
end

% -------------------- Plotting --------------------
if p.Results.plot

    if p.Results.newfig
        f = figure;
    end

    if p.Results.raincloud
        tiledlayout(1, nCond + 1, ...
            'TileSpacing', 'compact', ...
            'Padding', 'compact');
    else
        tiledlayout(1, nCond, ...
            'TileSpacing', 'compact', ...
            'Padding', 'compact');
    end

    % ----- Shared map limits -----
    maplimits = p.Results.maplimits;

    if isempty(maplimits)
        allVals = [];

        for c = 1:nCond
            allVals = [allVals; topoData{c}(:)];
        end

        allVals = allVals(isfinite(allVals));
        absMax = max(abs(allVals(:)));

        if isempty(absMax) || ~isfinite(absMax) || absMax == 0
            maplimits = 'absmax';
        else
            maplimits = [-absMax absMax];
        end
    end

    % ----- Topoplots -----
    for c = 1:nCond

        nexttile(c);

        thisTopoplotArgs = p.Results.topoplotArgs;

        if p.Results.markChannels && ~isempty(chanIdx)
            thisTopoplotArgs = [thisTopoplotArgs, ...
                {'emarker2', {chanIdx, char(p.Results.markMarker), p.Results.markColor, ...
                p.Results.markSize, p.Results.markLineWidth}}];
        end
        
        topoplot(topoData{c}, chanlocs, ...
            'maplimits', maplimits, ...
            thisTopoplotArgs{:});
        colormap(p.Results.cmap);

        if doFilter
            title(sprintf('%s | %g-%g Hz | %g-%g ms', ...
                names{c}, freqWindow(1), freqWindow(2), timeWindow(1), timeWindow(2)), ...
                'Interpreter', 'none');
        else
            title(sprintf('%s | %g-%g ms', ...
                names{c}, timeWindow(1), timeWindow(2)), ...
                'Interpreter', 'none');
        end

        colorbar;
    end

    % ----- Raincloud plot -----
    if p.Results.raincloud

        % IMPORTANT:
        % Do NOT call nexttile manually here.
        % bbar_raincloud handles the tile when given 'nexttile'.
        [~, Stats] = bbar_raincloud(rcData, ...
            'nexttile', nCond + 1, ...
            'xtickslabels', [names], ...
            'ylabel', sprintf('%s at selected channel(dB)', measure), ...
            p.Results.raincloudArgs{:});
    end

    set(gcf, 'Color', [1 1 1]);
end

end

% ========================================================================
% Helper: find channels
% ========================================================================
function chanIdx = find_channels(chan, chanlocs)

labels = string({chanlocs.labels});

if isnumeric(chan)
    chanIdx = chan;
    return
end

if ischar(chan) || isstring(chan)
    chan = cellstr(string(chan));
end

if iscell(chan)
    chanIdx = nan(1, numel(chan));

    for i = 1:numel(chan)
        hit = find(strcmpi(labels, string(chan{i})), 1);

        if isempty(hit)
            error('Channel %s was not found.', string(chan{i}));
        end

        chanIdx(i) = hit;
    end
else
    error('chan must be numeric, char, string, or cell array of channel labels.');
end

end