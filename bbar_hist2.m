
function h = bbar_hist2(X, Y, xbins, ybins)
    % 2D histogram computation
    xmin = min(X(:));
    xmax = max(X(:));
    ymin = min(Y(:));
    ymax = max(Y(:));
    
    % Create histogram
    h = zeros(xbins, ybins);
    
    % Bin edges
    xedges = linspace(xmin, xmax, xbins+1);
    yedges = linspace(ymin, ymax, ybins+1);
    
    % Compute 2D histogram
    for i = 1:numel(X)
        xbin = find(X(i) >= xedges(1:end-1) & X(i) < xedges(2:end));
        ybin = find(Y(i) >= yedges(1:end-1) & Y(i) < yedges(2:end));
        
        if ~isempty(xbin) && ~isempty(ybin)
            h(xbin, ybin) = h(xbin, ybin) + 1;
        end
    end
end