function [FullStats, erpData] = bbar_plotERP(data, varargin)
% Parser
p = inputParser;
addRequired(p, 'data', @iscell);
addRequired(p, 'chan', @isstring);
addOptional(p, 'plot', 1, @isnumeric);
addOptional(p, 'plotstd', 1, @isnumeric);
addOptional(p, 'grid', 1, @isnumeric);
addOptional(p, 'newfig', 1, @isnumeric);
addOptional(p, 'ylim',[], @isnumeric);
addOptional(p, 'xlabel', 'time (ms)', @ischar);
addOptional(p, 'ylabel', 'μV', @ischar);
addOptional(p, 'markedarea', NaN, @isnumeric);
addOptional(p, 'baseline', NaN, @isnumeric);
addOptional(p, 'myColormap', bbar_makecmap([51/255, 92/255, 103/255; 224/255, 159/255, 62/255; 84/255, 11/255, 14/255], size(data,2)));
addOptional(p, 'raincloud', [], @(x) isnumeric(x) && (isempty(x) || numel(x)==2));
addOptional(p, 'raincloudArgs', {}, @iscell);
addOptional(p, 'legend', 1, @isnumeric);
parse(p, data, varargin{:});

% -------------------- Extract time vector ONCE --------------------
% This is the authoritative time axis used for all indexing and plotting
times = data{1}(1).times;   % 1 x nTimepoints, in ms

% -------------------- Baseline handling --------------------
doBaseline = ~(isempty(p.Results.baseline) || (isnumeric(p.Results.baseline) && any(isnan(p.Results.baseline))));
if doBaseline
    bl = p.Results.baseline;
    if ~isnumeric(bl) || numel(bl) ~= 2
        error('baseline must be a 1x2 numeric vector in ms, e.g. [-200 0], or NaN/[] to disable.');
    end
    bl = sort(bl(:))';
    blIdx = find(times >= bl(1) & times <= bl(2));
    if isempty(blIdx)
        error('Baseline window [%g %g] ms does not overlap with times range [%g %g] ms.', ...
            bl(1), bl(2), times(1), times(end));
    end
end

% -------------------- Main code --------------------
numberOfDatasets = size(data,2);
uniqueSbs = unique(string({data{1}.subject}'));
erpData = cell(size(data));   % each cell: nSubs x nTimepoints

FullStats = table();

% Fetch dataset names from the EEG setname field
datasetNames = cell(1, numberOfDatasets);
for i1 = 1:numberOfDatasets
    try
        datasetNames{i1} = data{i1}(1).setname;
    catch
        datasetNames{i1} = sprintf('Condition %d', i1);
    end
end

for i1 = 1:numberOfDatasets
    k = 1;
    for i2 = uniqueSbs'
        chanIdx = find(strcmp(string({data{i1}(1).chanlocs.labels}'), p.Results.chan));
        idx = find(strcmp(string({data{1}.subject}), i2));
        concatenatedData = cat(3, data{i1}(idx).data);          % chan x nTimepoints x nTrials
        averagedData = mean(concatenatedData(chanIdx,:,:), 3);   % 1 x nTimepoints
        if doBaseline
            blMean = mean(averagedData(blIdx), 2);               % scalar, indexed by blIdx into nTimepoints
            averagedData = averagedData - blMean;
        end
        erpData{i1}(k,:) = averagedData;                        % nSubs x nTimepoints
        k = k + 1;
    end
end

% -------------------- ylim: include SEM extent --------------------
xtrems = [];
if isempty(p.Results.ylim)
    for i = 1:length(erpData)
        semCurve = std(erpData{i}, 1) / sqrt(size(erpData{i}, 1));  % 1 x nTimepoints
        upperEnv = max(mean(erpData{i}, 1) + semCurve);
        lowerEnv = min(mean(erpData{i}, 1) - semCurve);
        xtrems = [xtrems; max(abs([upperEnv, lowerEnv]))];
    end
    myYlim = [-max(xtrems) max(xtrems)];
else
    myYlim = p.Results.ylim;
end

% -------------------- Plotting --------------------
doRaincloud = ~isempty(p.Results.raincloud);

if p.Results.plot
    if p.Results.newfig
        f = figure;
    end

    if doRaincloud
        tl = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
        ax_erp = nexttile(1, [1 2]);
    else
        ax_erp = gca;
    end

    % ----- ERP axes -----
    legendHandles = gobjects(numberOfDatasets, 1);
    for i1 = 1:size(erpData,2)
    
        % Marked area: find time indices from times vector, then use actual ms values
        if ~isnan(p.Results.markedarea)
            maIdx = find(times >= p.Results.markedarea(1) & times <= p.Results.markedarea(2));
            if ~isempty(maIdx)
                x = [times(maIdx(1))  times(maIdx(1))  times(maIdx(end))  times(maIdx(end))];
                y = [myYlim(1) myYlim(2) myYlim(2) myYlim(1)];
                fill(x, y, 'k', 'FaceAlpha', 0.05, 'EdgeColor', 'none');
            else
                warning('bbar_plotERP: markedarea [%g %g] ms does not overlap with times [%g %g] ms. Skipping shading.', ...
                    p.Results.markedarea(1), p.Results.markedarea(2), times(1), times(end));
            end
        end
    
        % SEM envelope: x-axis is times vector
        if p.Results.plotstd
            curve1 = mean(erpData{i1},1) + (std(erpData{i1},1) / sqrt(size(erpData{i1},1)));
            curve2 = mean(erpData{i1},1) - (std(erpData{i1},1) / sqrt(size(erpData{i1},1)));
            x2 = [times, fliplr(times)];
            inBetween = [curve1, fliplr(curve2)];
            fill(x2, inBetween, p.Results.myColormap(i1,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end
    
        hold on
    
        % ERP line: x-axis is times vector
        legendHandles(i1) = plot(times,mean(erpData{i1},1), ...
            'Color', p.Results.myColormap(i1,:), 'LineWidth', 1.5);
    
        xline(0); yline(0);
        ylim(myYlim);
        xlim([times(1) times(end)]);       % full epoch range from times
    
        % xlabel anchored to last time point and bottom of ylim
        xl = xlabel(p.Results.xlabel);
        xl.Position(1) = times(end);
        xl.Position(2) = myYlim(1);
    
        ylabel(p.Results.ylabel);
        ax = gca;
        ax.XAxisLocation = "origin";
        ax.YAxisLocation = "origin";
        ax.Box = "off";
        if p.Results.grid; grid on; end
    end

    % Legend
    if p.Results.legend
        lgd = legend(legendHandles, datasetNames, 'Location', 'best');
        lgd.Box = 'off';
    end

    f.Color = [1 1 1];

    % ----- Raincloud axes -----
    if doRaincloud
        rcWin = sort(p.Results.raincloud);                      % [tStart tEnd] in ms
        rcIdx = find(times >= rcWin(1) & times <= rcWin(2));   % index into nTimepoints using times

        if isempty(rcIdx)
            warning('bbar_plotERP: raincloud window [%g %g] ms does not overlap with times [%g %g] ms. Skipping raincloud.', ...
                rcWin(1), rcWin(2), times(1), times(end));
        else
            % erpData{i1} is nSubs x nTimepoints → mean over rcIdx cols → nSubs x 1
            rcData = cell(1, numberOfDatasets);
            for i1 = 1:numberOfDatasets
                rcData{i1} = mean(erpData{i1}(:, rcIdx), 2);   % mean over time → nSubs x 1
            end

            ax_rc = nexttile(3);
            axes(ax_rc); %#ok<LAXES>
            [~, FullStats] = bbar_raincloud(rcData, ...
                'nexttile',     ax_rc, ...
                'cmap',         p.Results.myColormap, ...
                'ylabel',       'μV', ...
                'xlabel',       sprintf('%g – %g ms', rcWin(1), rcWin(2)), ...
                'xtickslabels', [{''},  datasetNames, {''}], ...
                p.Results.raincloudArgs{:});
        end
    end
end
end