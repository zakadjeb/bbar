% [idx,amp] = bbar_findpeaks(data, posneg)
%
% Find positive or negative peaks
%
%
% data                    - (Required) Data as vector, e.g. 1000 time-series.
%
% method                  - (Optional) String for method, e.g. 'diff' or
%                           'gradient'
%
% posneg                  - (Optional) Find positive or negative peaks
%                           (default: 'positive')
%
% Outputs:
%   idx                   - The index of the input signal
%   amp                   - The amplitude of the input signal
%
% See also:
%   bbar_raincloud, bbar_makecmap
%
%
% Zak Djebbara, Oct. 2023
%
function [idx,amp] = bbar_findpeaks(data,method,posneg)
if ~isnumeric(data) || (size(data,1) < 1 && size(data,2) < 1); error('Use a vector as input'); end

if nargin == 1; method = 'diff'; posneg = 'positive'; end
if nargin == 2; posneg = 'positive'; end

if size(data,1) > size(data,2); data = data'; end

if strcmp(method,'diff')
    data_dx = diff(data);
elseif strcmp(method,'gradient')
    data_dx = gradient(data);
end


if strcmp(posneg,'positive')
    % Identifying regular peaks
    idx = find(data_dx(2:end) < 0 & data_dx(1:end-1) > 0);

    % Making sure to get the initial point as a peak
    if data(2)-data(1) < 0
        idx = [0, idx];
    end
    % Making sure to get the last point as a peak
    if data(end)-data(end-1) > 0
        idx = [idx, length(data)-1];
    end

    idx = idx + 1;
    amp = data(idx);
elseif strcmp(posneg,'negative')
    % Identifying regular peaks
    idx = find(data_dx(1:end-1) < 0 & data_dx(2:end) > 0);

    % Making sure to get the initial point as a peak
    if data(2)-data(1) > 0
        idx = [0, idx];
    end
    % Making sure to get the last point as a peak
    if data(end)-data(end-1) < 0
        idx = [idx, length(data)-1];
    end
    
    idx = idx + 1;
    amp = data(idx);
elseif strcmp(posneg,'mean')
    idx = [];
    amp = mean(data);
end

if isempty(amp) && strcmp(posneg,'negative')
    midWay = round(length(data)/2);
    amp = data(midWay);
    % amp = min(data);
    idx = find(amp);
elseif isempty(amp) && strcmp(posneg,'positive')
    midWay = round(length(data)/2);
    amp = data(midWay);
    % amp = max(data);
    idx = find(amp);
end
end