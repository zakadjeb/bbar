function [nOpt, energyCum, s] = bbar_dmd_determine_modes(X, tau, nmax)
% bbar_dmd_determine_modes  Choose DMD rank using SVD energy only.
%
%   [nOpt, energyCum, s] = bbar_dmd_determine_modes(X, tau, nmax)
%
% INPUT
%   X     : [features x time] snapshot matrix (e.g., pixels x frames)
%   tau   : energy threshold in (0,1], e.g. 0.95 (default)
%   nmax  : optional cap on modes to consider (default: min(size(X))-1)
%
% OUTPUT
%   nOpt      : smallest n such that cumulative SVD energy >= tau
%   energyCum : cumulative energy curve (same length as sUsed)
%   s         : singular values used (vector)
%
% Notes:
%   Uses x1 = X(:,1:end-1), as in standard DMD.
%   Energy is computed from squared singular values.

if nargin < 2 || isempty(tau),  tau  = 0.95; end
if nargin < 3 || isempty(nmax), nmax = inf;  end

% DMD-style snapshot matrix
x1 = X(:,1:end-1);

% Econ SVD singular values
s = svd(x1,'econ');

% Cap number of modes considered
ncap = min([numel(s), nmax]);
s = s(1:ncap);

% Cumulative energy
energyCum = cumsum(s.^2) ./ sum(s.^2);

% Pick smallest n meeting threshold
nOpt = find(energyCum >= tau, 1, 'first');
if isempty(nOpt), nOpt = ncap; end
end
