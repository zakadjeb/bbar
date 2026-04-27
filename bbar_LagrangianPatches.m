function out = bbar_LagrangianPatches(videoPath, varargin)
% Extract Lagrangian (trajectory-aligned) patches and patch time series.
%
% OUTPUT:
% out.pos{t}   : Nx2 positions [x y] at frame t (1-indexed frames)
% out.signal   : NxT signal (mean intensity of each patch per frame)
% out.valid    : NxT logical, true if patch was fully inside frame
% out.params   : struct of parameters used
%
% Notes:
% - This does NOT stabilize the scene (good for retinal-motion focus).
% - Patches follow optical-flow advection.
%
% Zakaria-style: keep it modular and inspect intermediate products.

p = inputParser;
addRequired(p, 'videoPath', @(x)ischar(x) || isstring(x));

addParameter(p, 'SeedStride', 20, @(x)isnumeric(x) && isscalar(x) && x>=1);
addParameter(p, 'PatchRadius', 8, @(x)isnumeric(x) && isscalar(x) && x>=1);  % patch size = (2r+1)^2
addParameter(p, 'MaxFrames', inf, @(x)isnumeric(x) && isscalar(x) && x>=1);
addParameter(p, 'UseCLAHE', false, @(x)islogical(x) && isscalar(x));
addParameter(p, 'SignalType', 'mean', @(x)ischar(x) || isstring(x)); % 'mean'|'median'|'rmscontrast'
addParameter(p, 'FlowModel', 'farneback', @(x)ischar(x) || isstring(x)); % placeholder for extensions
parse(p, videoPath, varargin{:});
S = p.Results;

vr = VideoReader(videoPath);

% Read first frame
if ~hasFrame(vr)
    error('Video has no frames.');
end
I0 = readFrame(vr);
I0g = toGray(I0, S.UseCLAHE);

H = size(I0g,1);
W = size(I0g,2);

% Seed points on a grid
[xGrid, yGrid] = meshgrid( (1+S.PatchRadius):S.SeedStride:(W-S.PatchRadius), ...
                           (1+S.PatchRadius):S.SeedStride:(H-S.PatchRadius) );
x = xGrid(:);
y = yGrid(:);
N = numel(x);

% Prepare storage (we don't know T in advance if MaxFrames=inf)
pos = {};  % pos{t} = [x y] positions
pos{1} = [x y];

% For signal, we'll store in a growable array then trim
signal = nan(N, 1);
valid  = false(N, 1);

% Optical flow object
of = opticalFlowFarneback;

% Initialize flow with first frame
estimateFlow(of, I0g);

% Compute signal at t=1
[signal(:,1), valid(:,1)] = patchSignalAt(I0g, pos{1}, S.PatchRadius, S.SignalType);

t = 1;
while hasFrame(vr) && t < S.MaxFrames
    t = t + 1;
    I = readFrame(vr);
    Ig = toGray(I, S.UseCLAHE);

    flow = estimateFlow(of, Ig); % flow between previous and current frame internally

    % Advect points using current flow field sampled at previous positions
    prevPos = pos{t-1};
    [u, v] = sampleFlow(flow, prevPos(:,1), prevPos(:,2), W, H);

    newX = prevPos(:,1) + u;
    newY = prevPos(:,2) + v;

    pos{t} = [newX newY];

    % Extract patch signals at new positions
    [signal(:,t), valid(:,t)] = patchSignalAt(Ig, pos{t}, S.PatchRadius, S.SignalType);
end

out.pos    = pos;
out.signal = signal;
out.valid  = valid;
out.params = S;
end

% ----------------- helpers -----------------

function Ig = toGray(I, useCLAHE)
if size(I,3) == 3
    Ig = rgb2gray(I);
else
    Ig = I;
end
Ig = im2single(Ig);
if useCLAHE
    Ig = adapthisteq(Ig);
end
end

function [u, v] = sampleFlow(flow, x, y, W, H)
% flow.Vx, flow.Vy are HxW matrices
% Sample with bilinear interpolation; clamp coords.
xq = min(max(x, 1), W);
yq = min(max(y, 1), H);

u = interp2(flow.Vx, xq, yq, 'linear', 0);
v = interp2(flow.Vy, xq, yq, 'linear', 0);
end

function [sig, isValid] = patchSignalAt(Ig, posXY, r, signalType)
H = size(Ig,1); W = size(Ig,2);
x = posXY(:,1); y = posXY(:,2);

% Patch bounds
x1 = round(x - r); x2 = round(x + r);
y1 = round(y - r); y2 = round(y + r);

isValid = x1>=1 & y1>=1 & x2<=W & y2<=H;
sig = nan(size(x));

% Only compute for valid patches
idx = find(isValid);
if isempty(idx), return; end

ps = 2*r + 1;

switch lower(string(signalType))
    case "mean"
        for k = 1:numel(idx)
            i = idx(k);
            P = Ig(y1(i):y2(i), x1(i):x2(i));
            sig(i) = mean(P(:));
        end

    case "median"
        for k = 1:numel(idx)
            i = idx(k);
            P = Ig(y1(i):y2(i), x1(i):x2(i));
            sig(i) = median(P(:));
        end

    case "rmscontrast"
        % RMS contrast = std / mean (or just std), choose what you want.
        for k = 1:numel(idx)
            i = idx(k);
            P = Ig(y1(i):y2(i), x1(i):x2(i));
            mu = mean(P(:));
            sd = std(P(:));
            sig(i) = sd / max(mu, 1e-6);
        end

    otherwise
        error('Unknown SignalType: %s', signalType);
end
end