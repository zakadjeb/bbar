% bbar_makecmap()    - Creates a linearly related colormap based on a (N x 3) matrix.
%                         
%
% Usage:
%   >> [cmap] =   bbar_makecmap( colormatrix , numcolors)
%   >> [cmap] =   bbar_makecmap( [0,0,0; 1,1,1], 1000);
%   >> size(cmap) = 1000 x 3
%
% hint: use hex2rgb('ccccc')/255
%
% Inputs:
%   colormatrix             - [N x 3] color matrix, e.g. 
%                             [255 255 255; 150 150 150; 100 100 100]
%   numcolors               - number of colors to produce.
%   loop                    - whether or not to loop instead of creating a
%                             linear colormap (DEFAULT: 0)
%
% Outputs:
%   cmap                    - interpolated color matrix
%
%
% Author: Zak Djebbara, 2022
function [cmap] = bbar_makecmap(incmap, numcmap, loop)
if nargin == 0; help bbar_makecmap; end
if nargin <= 2; loop = 0; end
if size(incmap,2) < 3 || size(incmap,1) < 2; error('Input does not match criteria. Please make sure the input is at least 2 x 3 matrix'); end
if size(incmap,1) == numcmap; cmap = incmap; return; end
if ~isnumeric(numcmap); error('Please make sure the numbcmap variable is a number'); end

totalnum = numcmap-size(incmap, 1);
if totalnum<0
    warning('There are more colors inserted than wished to generate.');
    cmap = incmap(1:numcmap,:);
    return;
end

cmap = cell(numcmap,1);
jumpIdx = ceil(numcmap/(size(incmap,1)-1));
serie = 1:jumpIdx-1:numcmap;
if rem(numcmap,size(incmap,1)) ~= 0 && numcmap-serie(end) == 1; serie(end+1) = size(cmap,1); end
if length(serie) > size(incmap,1); serie(end-1) = []; end
if length(serie) < size(incmap,1); serie(end+1) = numcmap; end
% if serie(end) ~= numcmap; serie(end) = numcmap; end
while numcmap-serie(end) ~= 0
    if numcmap-serie(end) < jumpIdx
        serie(end) = serie(end) + 1;
    elseif numcmap-serie(end) > jumpIdx
        serie(end) = serie(end) - 1;
    end
end
for i = 1:length(serie)
    cmap{serie(i)} = incmap(i,:);
end



% Make so that it loops around the input numbers
if loop
    cmap = repmat(incmap, numcmap, 1);
    cmap = cmap(1:numcmap,:);
    return;
end

% Creating the new colors
tweennum = ((serie(2:end)-serie(1:end-1))-1)';
amountoftimes = size(tweennum,1);
newColors = cell(1,amountoftimes);
for i1 = 1:amountoftimes
    newr = []; newg = []; newb = [];
    newr = interp1(incmap(i1:i1+1,1)', linspace(1,2, tweennum(i1)+2));
    newr = newr(2:end-1)';
    newg = interp1(incmap(i1:i1+1,2)', linspace(1,2, tweennum(i1)+2));
    newg  = newg (2:end-1)';
    newb = interp1(incmap(i1:i1+1,3)', linspace(1,2, tweennum(i1)+2));
    newb = newb(2:end-1)';
    newColors{i1} = [newr, newg, newb];
    
end

% Making sure the newColors fit the identified idx
newColors = unpackCell(newColors);
idx = findemptycell(cmap);
if size(newColors,1) ~= length(idx); error('something went wrong.'); end

% Inserting the newColors according to the identified idx
for i = 1:size(idx,1)
    cmap{idx(i)} = newColors(i,:);
end

% Preparing output
cmap = cell2mat(cmap);

end

function idx = findemptycell(x)
idx = [];
for i = 1:length(x)
    if isempty(x{i})
        idx = [idx; i];
    end
end
end

function output = unpackCell(x)
output=[];
for i = 1:size(x,2)
    output = [output; x{i}];
end
end