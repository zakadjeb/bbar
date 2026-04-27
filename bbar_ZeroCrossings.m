function idx = bbar_ZeroCrossings(x)
% FINDZEROCROSSINGS  return indices where signal flips sign
%
% Example:
%   x = [2 1 -1 -3 4];
%   idx = findZeroCrossings(x)   % -> [2 4]

    x = x(:);                 % force column
    s = sign(x);              % +1 0 -1
    s(s==0) = nan;            % remove exact zeros, so they don't count as "no sign"
    ds = diff(s);
    % a sign flip is where the diff goes nonzero
    idx = find(~isnan(ds) & abs(ds) == 2); 
end