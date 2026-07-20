function [FullStats, psdData, freqs, bandData] = bbar_psdplot(data, chan, varargin)
% BBAR_PSDPLOT Plot condition-wise PSDs and optionally run raincloud stats.
%
%   [FullStats, psdData, freqs, bandData] = bbar_psdplot(data, chan, ...)
%
% INPUT
%   data : cell array where each cell contains an EEGLAB-like EEG struct array.
%          Each element in data{i} is treated as one subject/observation.
%          EEG.data must be channels x timepoints x trials, or channels x timepoints.
%   chan : channel label, e.g. "Cz" or 'Cz'.
%
% NAME-VALUE INPUTS
%   'plot'           : 1/0, whether to plot. Default = 1.
%   'plotstd'        : 1/0, whether to plot SEM envelope. Default = 1.
%   'grid'           : 1/0, show grid. Default = 1.
%   'newfig'         : 1/0, create new figure. Default = 1.
%   'ylim'           : y-axis limits. Default = [].
%   'xlim'           : x-axis frequency limits in Hz. Default = [].
%   'xlabel'         : x-axis label. Default = 'Frequency (Hz)'.
%   'ylabel'         : y-axis label. Default = 'PSD (dB/Hz)'.
%   'timeWindow'     : optional time window in ms before PSD, e.g. [300 800].
%                      Requires EEG.times. Default = [].
%   'freqWindow'     : frequency window in Hz for shading + raincloud/stat extraction,
%                      e.g. [8 12]. Default = [].
%   'raincloud'      : backwards-compatible alias for 'freqWindow'. If supplied, it
%                      overrides 'freqWindow'. Default = [].
%   'bandMeasure'    : value extracted per subject from freqWindow:
%                      'mean', 'median', 'sum', 'trapz', 'max', 'peak', or 'min'.
%                      Default = 'mean'.
%   'psdUnits'       : 'db' or 'linear'. Default = 'db'.
%   'welchWindow'    : pwelch window length in samples, or a vector window. Default = [].
%   'noverlap'       : pwelch overlap in samples. Default = [].
%   'nfft'           : pwelch nfft. Default = [].
%   'detrend'        : 1/0, detrend each trial before PSD. Default = 1.
%   'myColormap'     : nConditions x 3 RGB colormap. Default = bbar_makecmap(...).
%   'raincloudArgs'  : cell array passed directly to bbar_raincloud for statistics.
%   'legend'         : 1/0, show legend. Default = 1.
%
% OUTPUT
%   FullStats : table returned by bbar_raincloud, if freqWindow/raincloud is supplied.
%   psdData   : cell array. psdData{i} is nSubjects x nFreqs.
%   freqs     : PSD frequency vector in Hz.
%   bandData  : cell array. bandData{i} is one extracted band value per subject.
%
% EXAMPLE
%   [Stats, psdData, freqs, alphaPower] = bbar_psdplot(data, "Cz", ...
%       'freqWindow', [8 12], ...
%       'xlim', [1 40], ...
%       'bandMeasure', 'mean', ...
%       'raincloudArgs', {'test', 'friedman'});

% -------------------- Parser --------------------
p = inputParser;
p.FunctionName = 'bbar_psdplot';

addRequired(p, 'data', @iscell);
addRequired(p, 'chan', @(x) ischar(x) || isstring(x));

addParameter(p, 'plot', 1, @isnumeric);
addParameter(p, 'plotstd', 1, @isnumeric);
addParameter(p, 'grid', 1, @isnumeric);
addParameter(p, 'newfig', 1, @isnumeric);
addParameter(p, 'ylim', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'xlim', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 2));
addParameter(p, 'xlabel', 'Frequency (Hz)', @(x) ischar(x) || isstring(x));
addParameter(p, 'ylabel', 'PSD (dB/Hz)', @(x) ischar(x) || isstring(x));
addParameter(p, 'timeWindow', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 2));
addParameter(p, 'freqWindow', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 2));
addParameter(p, 'raincloud', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 2));
addParameter(p, 'bandMeasure', 'mean', @(x) ischar(x) || isstring(x));
addParameter(p, 'psdUnits', 'db', @(x) ischar(x) || isstring(x));
addParameter(p, 'welchWindow', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'noverlap', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'nfft', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'detrend', 1, @isnumeric);
addParameter(p, 'myColormap', bbar_makecmap([51/255, 92/255, 103/255; 224/255, 159/255, 62/255; 84/255, 11/255, 14/255], numel(data)), @isnumeric);
addParameter(p, 'raincloudArgs', {}, @iscell);
addParameter(p, 'legend', 1, @isnumeric);

parse(p, data, chan, varargin{:});

bandMeasure = lower(char(p.Results.bandMeasure));
validBandMeasures = {'mean','median','sum','trapz','max','peak','min'};
if ~ismember(bandMeasure, validBandMeasures)
    error('bandMeasure must be ''mean'', ''median'', ''sum'', ''trapz'', ''max'', ''peak'', or ''min''.');
end

psdUnits = lower(char(p.Results.psdUnits));
validPsdUnits = {'db','linear'};
if ~ismember(psdUnits, validPsdUnits)
    error('psdUnits must be ''db'' or ''linear''.');
end

% raincloud is kept as a backwards-compatible alias for the analysis window.
% If both are supplied, raincloud takes precedence.
if ~isempty(p.Results.raincloud)
    statWin = sort(p.Results.raincloud(:))';
elseif ~isempty(p.Results.freqWindow)
    statWin = sort(p.Results.freqWindow(:))';
else
    statWin = [];
end

doBandStats = ~isempty(statWin);

doTimeWindow = ~isempty(p.Results.timeWindow);
if doTimeWindow
    tw = sort(p.Results.timeWindow(:))';
end

% -------------------- Main PSD extraction --------------------
numberOfDatasets = numel(data);
psdData = cell(size(data));
bandData = cell(size(data));
FullStats = table();
freqs = [];

% Fetch dataset names from the EEG setname field.
datasetNames = cell(1, numberOfDatasets);
for i1 = 1:numberOfDatasets
    try
        datasetNames{i1} = data{i1}(1).setname;
    catch
        datasetNames{i1} = sprintf('Condition %d', i1);
    end
end

for i1 = 1:numberOfDatasets

    nRows = numel(data{i1});
    if nRows == 0
        warning('bbar_psdplot: data{%d} is empty. Skipping this dataset.', i1);
        psdData{i1} = [];
        continue
    end

    if ~isfield(data{i1}(1), 'chanlocs') || ~isfield(data{i1}(1).chanlocs, 'labels')
        error('data{%d}(1).chanlocs.labels was not found.', i1);
    end

    chanIdx = find(strcmp(string({data{i1}(1).chanlocs.labels}'), string(p.Results.chan)));
    if isempty(chanIdx)
        error('Channel %s was not found in data{%d}.', string(p.Results.chan), i1);
    elseif numel(chanIdx) > 1
        warning('bbar_psdplot: Multiple channels named %s found in data{%d}. Using the first match.', string(p.Results.chan), i1);
        chanIdx = chanIdx(1);
    end

    for iRow = 1:nRows

        if ~isfield(data{i1}(iRow), 'srate') || isempty(data{i1}(iRow).srate)
            error('data{%d}(%d).srate is missing. PSD analysis requires the sampling rate.', i1, iRow);
        end
        fs = double(data{i1}(iRow).srate);

        thisData = data{i1}(iRow).data;              % chan x nTimepoints x nTrials, or chan x nTimepoints
        thisChan = squeeze(thisData(chanIdx,:,:));  % nTimepoints x nTrials, usually

        % Force orientation to nTimepoints x nTrials.
        if isvector(thisChan)
            thisChan = thisChan(:);
        else
            % After squeeze, EEGLAB epoch data should already be time x trials.
            % If it is trials x time, transpose it when EEG.pnts gives the clue.
            if isfield(data{i1}(iRow), 'pnts') && size(thisChan, 2) == data{i1}(iRow).pnts && size(thisChan, 1) ~= data{i1}(iRow).pnts
                thisChan = thisChan';
            end
        end

        % Optional epoch/time cropping before PSD.
        if doTimeWindow
            if ~isfield(data{i1}(iRow), 'times') || isempty(data{i1}(iRow).times)
                error('timeWindow was supplied, but data{%d}(%d).times is missing.', i1, iRow);
            end
            times = double(data{i1}(iRow).times(:));
            tIdx = find(times >= tw(1) & times <= tw(2));
            if isempty(tIdx)
                error('timeWindow [%g %g] ms does not overlap with data{%d}(%d).times [%g %g] ms.', ...
                    tw(1), tw(2), i1, iRow, times(1), times(end));
            end
            thisChan = thisChan(tIdx, :);
        end

        % Remove trials that are completely NaN.
        goodTrials = ~all(isnan(thisChan), 1);
        thisChan = thisChan(:, goodTrials);
        if isempty(thisChan)
            warning('bbar_psdplot: data{%d}(%d) has no usable trials after NaN removal.', i1, iRow);
            continue
        end

        if p.Results.detrend
            thisChan = detrend(thisChan);
        end

        nSamples = size(thisChan, 1);

        % Welch parameters.
        if isempty(p.Results.welchWindow)
            winLength = min(nSamples, max(16, floor(nSamples / 2)));
            win = hamming(winLength);
        elseif isscalar(p.Results.welchWindow)
            winLength = min(nSamples, round(p.Results.welchWindow));
            win = hamming(winLength);
        else
            win = p.Results.welchWindow(:);
            winLength = numel(win);
            if winLength > nSamples
                error('welchWindow is longer than the available signal in data{%d}(%d).', i1, iRow);
            end
        end

        if isempty(p.Results.noverlap)
            noverlap = floor(0.5 * winLength);
        else
            noverlap = p.Results.noverlap;
        end

        if noverlap >= winLength
            error('noverlap must be smaller than the Welch window length.');
        end

        if isempty(p.Results.nfft)
            nfft = max(256, 2^nextpow2(winLength));
        else
            nfft = p.Results.nfft;
        end

        % Compute one PSD per trial, then average trials within subject/observation.
        [Pxx, f] = pwelch(thisChan, win, noverlap, nfft, fs);  % nFreqs x nTrials
        subjPsd = mean(Pxx, 2, 'omitnan')';                   % 1 x nFreqs

        if strcmp(psdUnits, 'db')
            subjPsd = 10 * log10(subjPsd);
        end

        % Use the first frequency vector as authoritative. If later data produce
        % another vector, interpolate onto the first one.
        if isempty(freqs)
            freqs = f(:)';
        elseif numel(f) ~= numel(freqs) || any(abs(f(:)' - freqs) > 1e-10)
            subjPsd = interp1(f(:)', subjPsd, freqs, 'linear', NaN);
        end

        psdData{i1}(iRow, :) = subjPsd;
    end
end

% -------------------- Extract band values for statistics --------------------
if doBandStats
    if isempty(freqs)
        error('No PSD frequency vector was produced. Check that data contain usable observations.');
    end

    fIdx = find(freqs >= statWin(1) & freqs <= statWin(2));
    if isempty(fIdx)
        warning('bbar_psdplot: freqWindow/raincloud [%g %g] Hz does not overlap with PSD frequencies [%g %g] Hz. Skipping band extraction.', ...
            statWin(1), statWin(2), freqs(1), freqs(end));
        doBandStats = false;
    else
        for i1 = 1:numberOfDatasets
            if isempty(psdData{i1})
                bandData{i1} = [];
                continue
            end

            thisBand = psdData{i1}(:, fIdx);  % nSubjects x nBandFreqs

            switch bandMeasure
                case 'mean'
                    bandData{i1} = mean(thisBand, 2, 'omitnan');
                case 'median'
                    bandData{i1} = median(thisBand, 2, 'omitnan');
                case {'sum','trapz'}
                    bandData{i1} = trapz(freqs(fIdx), thisBand, 2);
                case {'max','peak'}
                    bandData{i1} = max(thisBand, [], 2, 'omitnan');
                case 'min'
                    bandData{i1} = min(thisBand, [], 2, 'omitnan');
            end
        end
    end
end

% -------------------- ylim: include SEM extent --------------------
if isempty(p.Results.ylim)
    xtrems = [];
    for i1 = 1:numberOfDatasets
        if isempty(psdData{i1})
            continue
        end
        meanCurve = mean(psdData{i1}, 1, 'omitnan');
        nAtFreq = sum(~isnan(psdData{i1}), 1);
        semCurve = std(psdData{i1}, 0, 1, 'omitnan') ./ sqrt(nAtFreq);
        upperEnv = max(meanCurve + semCurve, [], 'omitnan');
        lowerEnv = min(meanCurve - semCurve, [], 'omitnan');
        xtrems = [xtrems; upperEnv; lowerEnv]; %#ok<AGROW>
    end

    if isempty(xtrems) || all(isnan(xtrems))
        myYlim = [-1 1];
    else
        pad = 0.05 * range(xtrems);
        if pad == 0 || isnan(pad)
            pad = 1;
        end
        myYlim = [min(xtrems) - pad, max(xtrems) + pad];
    end
else
    myYlim = p.Results.ylim;
end

if isempty(p.Results.xlim)
    myXlim = [freqs(1), freqs(end)];
else
    myXlim = sort(p.Results.xlim(:))';
end

% -------------------- Plotting --------------------
if p.Results.plot
    if p.Results.newfig
        figure;
    end

    if doBandStats
        tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
        ax_psd = nexttile(1, [1 2]);
    else
        ax_psd = gca;
    end

    axes(ax_psd); %#ok<LAXES>
    hold(ax_psd, 'on');

    % Mark frequency band once, behind all curves.
    if doBandStats
        x = [statWin(1), statWin(1), statWin(2), statWin(2)];
        y = [myYlim(1), myYlim(2), myYlim(2), myYlim(1)];
        fill(ax_psd, x, y, 'k', 'FaceAlpha', 0.05, 'EdgeColor', 'none');
    end

    legendHandles = gobjects(numberOfDatasets, 1);

    for i1 = 1:numberOfDatasets
        if isempty(psdData{i1})
            continue
        end

        meanCurve = mean(psdData{i1}, 1, 'omitnan');

        if p.Results.plotstd
            nAtFreq = sum(~isnan(psdData{i1}), 1);
            semCurve = std(psdData{i1}, 0, 1, 'omitnan') ./ sqrt(nAtFreq);
            curve1 = meanCurve + semCurve;
            curve2 = meanCurve - semCurve;
            x2 = [freqs, fliplr(freqs)];
            inBetween = [curve1, fliplr(curve2)];
            fill(ax_psd, x2, inBetween, p.Results.myColormap(i1,:), ...
                'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end

        legendHandles(i1) = plot(ax_psd, freqs, meanCurve, ...
            'Color', p.Results.myColormap(i1,:), 'LineWidth', 1.5);
    end

    xlabel(ax_psd, p.Results.xlabel);
    ylabel(ax_psd, p.Results.ylabel);
    xlim(ax_psd, myXlim);
    ylim(ax_psd, myYlim);
    ax_psd.Box = 'off';
    if p.Results.grid
        grid(ax_psd, 'on');
    end

    if p.Results.legend
        validHandles = isgraphics(legendHandles);
        lgd = legend(ax_psd, legendHandles(validHandles), datasetNames(validHandles), 'Location', 'best');
        lgd.Box = 'off';
    end

    set(gcf, 'Color', [1 1 1]);

    % ----- Raincloud/stat axes -----
    if doBandStats
        ax_rc = nexttile(3);
        axes(ax_rc); %#ok<LAXES>

        [~, FullStats] = bbar_raincloud(bandData, ...
            'nexttile',     ax_rc, ...
            'cmap',         p.Results.myColormap, ...
            'ylabel',       p.Results.ylabel, ...
            'xlabel',       sprintf('%s PSD: %g–%g Hz', bandMeasure, statWin(1), statWin(2)), ...
            'xtickslabels', [{''}, datasetNames, {''}], ...
            p.Results.raincloudArgs{:});
    end
end
end
