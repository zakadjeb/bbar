function y = bbar_subdivide(x, param, method)
% bbar_subdivide - Downsample a vector by blockwise averaging.
%
% This function subdivides a signal into windows and returns the mean of
% each window. All data are always included. If the final window would be
% shorter than the intended window size, it is merged into the previous
% window.
%
% INPUTS:
%   x     - input time series vector
%   param - either:
%             * window size        if method = 'windowSize'
%             * number of outputs  if method = 'nSamples'
%   mode  - 'windowSize' or 'nSamples'
%
% OUTPUT:
%   y     - downsampled vector of window means
%
% EXAMPLES:
%   y = bbar_subdivide(x, 100, 'windowSize');
%   y = bbar_subdivide(x, 20,  'nSamples');

    x = x(:);
    L = length(x);

    if L == 0
        error('Input x is empty.');
    end

    if ~isscalar(param) || param <= 0 || param ~= round(param)
        error('param must be a positive integer scalar.');
    end

    switch lower(method)
        case 'windowsize'
            winSize = param;

            if winSize >= L
                y = mean(x);
                return;
            end

            starts = 1:winSize:L;
            ends = [starts(2:end)-1, L];

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
                y = x;
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

    y = zeros(numel(starts),1);
    for i = 1:numel(starts)
        y(i) = mean(x(starts(i):ends(i)));
    end
end