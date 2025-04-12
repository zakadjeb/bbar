% bbar_barplot()    - Creates a bar plot.
%                         
%
% Usage:
%   >> [cmap, tmpPval, tmpStats] = bbar_barplot(data,varargin)
%
% Inputs:
%   data                    - (Required) Data in cells:
%                             data{1} = 1x2 cells <- Condition 1 {Hits, Misses}
%                             data{2} = 1x2 cells <- Condition 2 {Hits, Misses}
%   cmap                    - Determines the colormap of the plots (default: bemobil_makecmap())
%   nexttile                - In case of using 'tiledlayout', set this to
%                             true/1 so that it plots it as 'nexttile'
%   xticks                  - Set the ticks on xline (default: based on data length)
%   xtickslabels            - Set the labels, eg. {'', 'Pre', 'Post'}
%   xlim                    - Set the extremes of the x axis.
%   stats                   - Compute the statistics, 'chi2' (default: [])
%   plotLsline              - Plots the least-square line to show a mean
%                             (default: 1)
%   multcmp                 - If there's more than one comparison, it will
%                             automatically correct for multiple
%                             comparison. Determine the method here:
%                             'holm', 'hochberg', 'hommel', 
%                             'bonferroni', 'BH', 'BY', 'fdr'(default), 'sidak' or 'none'.
%   comparison              - Comparison matrix, eg. [1 2; 2 3; 1 3] which compares
%                             condition 1 with 2, and 2 with 3, and 1 with 3 (default: [])
%   StatsYPos               - The location of the statistics in the plot.
%                             Must include two values. First one determines
%                             the y-position. Second one determines the 
%                             distance between the conditions (if multiple), e.g. [3 0.1].
%   Alpha                   - Alpha level for statistics (default: .05)
%   ylabel                  - The unit of data (default: 'magnitude')
%   legend                  - Add legend to plot, eg. {'Hit','Miss'} (default: [])
%   legendBox               - Determine whether to plot or not plot the box
%                             of the legend. (default: 1)
%   legendLoc               - Determine where to plot the box of the legend.
%                             (default: 'best')
%   title                   - Add title to plot (default: '')
%   grid                    - Add grid to plot (default: 'on')
%   printP                  - Determine if exact p-values should be printed
%                             with significance stars (default: 0)
%
%
% Outputs:
%   cmap                    - The color map used in the plot
%   tmpVal                  - The p-values, if computed
%   tmpStats                - The chi2 table, if computed
%
% See also:
%   bemobil_plot_eyetracking, bemobil_plot_heatmap, bemobil_plot_patterns
%
% Author: Zak Djebbara, 2022
%
%
function [b, tmpStats] = bbar_barplot(data,varargin)

for i = 1:length(data)
    defaultXTicksLabels{1,i} = ['Condition ' num2str(i)];
end
defaultXTicks = 0:length(data)+1;

p = inputParser;
addRequired(p, 'data', @iscell);
addOptional(p, 'alpha', .05, @isnumeric);
addOptional(p, 'cmap', [], @isnumeric);
addOptional(p, 'nexttile', 0);
addOptional(p, 'xticks', defaultXTicks, @isnumeric);
addOptional(p, 'xtickslabels', defaultXTicksLabels, @iscell);
addOptional(p, 'xlim', [.5 length(data)+.75], @isnumeric);
addOptional(p, 'ylim', [], @isnumeric);
addOptional(p, 'writeRatio', 0, @isnumeric);
addOptional(p, 'writeBarNames', {}, @iscell);
addOptional(p, 'plotLsline', 1, @isnumeric);
addOptional(p, 'stats', '', @ischar);
addOptional(p, 'multcmp','fdr', @ischar);
addOptional(p, 'comparison', [], @isnumeric);
addOptional(p, 'StatsYPos',[],@(x) (sum(size(x)) == 3) && isnumeric(x));
addOptional(p, 'ylabel', 'ocurrences', @ischar);
addOptional(p, 'legend', [], @iscell);
addOptional(p, 'legendLoc', 'best', @isstring);
addOptional(p, 'legendBox', [], @isnumeric);
addOptional(p, 'title', '', @ischar);
addOptional(p, 'grid', 'on');
addOptional(p, 'printP',0,@isnumeric);
parse(p,data,varargin{:});

% Setting up initial parameters 
% Set the colormap
if isempty(p.Results.cmap)
    mycmap = bemobil_makecmap([51/255, 92/255, 103/255; 224/255, 159/255, 62/255; 84/255, 11/255, 14/255], length(data{1}));
else
    mycmap = p.Results.cmap;
end
% Allow using tiledlayout
if p.Results.nexttile == 1
    nexttile;
end
% Ensure size of data matches xticklabels
if size(p.Results.xtickslabels,2) ~= size(data,2); error('XTicks labels do not match the size of data.'); end

%---------------------------
% Main plotting
nCol = size(data{1},2);
% figure;
b = bar(reshape(cell2mat(data),nCol,[])','FaceColor','flat');
set(gca,'XTickLabel', p.Results.xtickslabels);
for i1 = 1:length(b)
    for i2 = 1:size(b(i1).CData,1)
        b(i1).CData(i2,:) = mycmap(i1,:);
    end
end

if p.Results.plotLsline
    hold on
    for i1 = 1:nCol
        % plot(b(i1).XEndPoints,b(i1).YEndPoints,'Color',cmap(i1,:),'LineWidth',1.5);
        scatter(b(i1).XEndPoints,b(i1).YEndPoints,"filled",'.');
    end
    h = lsline;
    for i1 = 1:length(h)
        % h(i1).Color = mycmap(i1,:);
        h(i1).Color = 'k';
        h(i1).LineStyle = ':';
        h(i1).LineWidth = 1;
    end
end
% % Computing the position for the p-value text, if computed
% maxVal = [];
% for i = 1:nCol
%     maxVal = [maxVal; b(i).YEndPoints];
% end
% maxVal = max(max(maxVal));
% for i = 1:size(data,2)
%     distanceSignficancePlot(i) = (maxVal/10)*i;
% end


if isempty(p.Results.StatsYPos)
    distanceSignficancePlot = [];
    ax1 = gca;
    for i = 1:size(p.Results.comparison,1)
        distanceSignficancePlot = [distanceSignficancePlot; ax1.YTick(end)+(ax1.YTick(end)-ax1.YTick(end-1))*i];
    end
else
    distanceSignficancePlot = [];
    ax1 = gca;
    for i = 1:size(p.Results.comparison,1)
        if i == 1; distanceSignficancePlot = [distanceSignficancePlot; p.Results.StatsYPos(1)]; end
        if i ~= 1; distanceSignficancePlot = [distanceSignficancePlot; distanceSignficancePlot(end)+p.Results.StatsYPos(2)]; end
    end
end

% distanceSignficancePlot = [];
% ax1 = gca;
% for i = 1:size(p.Results.comparison,1)
%     distanceSignficancePlot = [distanceSignficancePlot; ax1.YTick(end)+(ax1.YTick(end)-ax1.YTick(end-1))*i];
% end

tmpStats{6} = [];
pvalStr = cell(size(p.Results.comparison,1),1);
condNames = p.Results.xtickslabels;

%---------------------------
% Computing the chi-squared
tmpPval{1} = []; tmpStats{1} = [];
if strcmp(p.Results.stats, 'chi2')
    for i = 1:size(p.Results.comparison,1)
        % Computing Chi2
        [pVal, type, chi2stat, df] = computeChi2([data{p.Results.comparison(i,1)};data{p.Results.comparison(i,2)}],p.Results.alpha);
        pVal = round(pVal,4,'decimal');
        tmpStats{i,1} = type;
        tmpStats{i,2} = [tmpStats{i,2}; condNames(p.Results.comparison(i,1)) + " vs " + condNames(p.Results.comparison(i,2))];
        if strcmp(type, 'Chi-squared test'); tmpStats{i,3} = [tmpStats{i,3}; sum(sum(chi2stat))]; tci = 'NA'; end
        if strcmp(type, 'Fishers exact test'); tmpStats{i,3} = [tmpStats{i,3}; chi2stat.OddsRatio]; tci = chi2stat.ConfidenceInterval; end
        tmpStats{i,4} = [tmpStats{i,4}; [round(mean(data{p.Results.comparison(i,1)}),2) round(mean(data{p.Results.comparison(i,1)}),2)]];
        tmpStats{i,5} = [tmpStats{i,5}; pVal];

        % Calculate Cramér's V
        if strcmp(type, 'Chi-squared test')
            numRows = size([data{p.Results.comparison(i,1)};data{p.Results.comparison(i,2)}], 1);
            numCols = size([data{p.Results.comparison(i,1)};data{p.Results.comparison(i,2)}], 2);
            N = sum(sum([data{p.Results.comparison(i,1)};data{p.Results.comparison(i,2)}]))/numRows;
            minDim = min(numRows - 1, numCols - 1);  % To avoid division by zero
            V = sqrt(sum(sum(chi2stat)) / N * minDim);
            tmpStats{i,6} = [tmpStats{i,6}; V];
        else
            tmpStats{i,6} = [tmpStats{i,6}; NaN];
        end

        % [pVal, chi2stat] = computeChi2([data{p.Results.comparison(i,1)};data{p.Results.comparison(i,2)}],p.Results.alpha);
        % tmpStats{1} = [tmpStats{1}; sum(sum(chi2stat))];
        % tmpPval{1} = [tmpPval{1}; pVal];
        % disp(['p-val: ' num2str(tmpPval)]);
    end
    
    tmpStats = cell2table(tmpStats);
    tmpStats = renamevars(tmpStats,["tmpStats1","tmpStats2","tmpStats3","tmpStats4","tmpStats5","tmpStats6"], ...
    ["Stats test","Comparison","Chi/Odds ratio","Means","P-value","Cramer's V"]);

    % Plotting p-values
    % Correcting for multiple comparison if comparions are more than one.
    if strcmp(p.Results.multcmp,'none')
        adj_p = tmpStats.("P-value");
    else
        adj_p = pval_adjust(tmpStats.("P-value"),p.Results.multcmp);
        if size(p.Results.comparison,1) > 1
            warning(['p-values are corrected for multiple comparison using ' p.Results.multcmp '.']);
        end
    end
    tmpStats.("P-value") = adj_p;
    for i = 1:size(p.Results.comparison,1)
        if tmpStats.("P-value")(i) <= 0.1 && tmpStats.("P-value")(i) > 0.05
            pvalStr{i} = ['' pvalStr{i} ' .'];
        elseif tmpStats.("P-value")(i) <= 0.05 && tmpStats.("P-value")(i) > 0.01
            pvalStr{i} = [pvalStr{i} ' *'];
        elseif tmpStats.("P-value")(i) <= 0.01 && tmpStats.("P-value")(i) > 0.001
            pvalStr{i} = [pvalStr{i} ' **'];
        elseif tmpStats.("P-value")(i) <= 0.001 && tmpStats.("P-value")(i) > 0.0001
            pvalStr{i} = [pvalStr{i} ' ***'];
        elseif tmpStats.("P-value")(i) <= 0.0001
            pvalStr{i} = [pvalStr{i} ' ****'];
        else
            pvalStr{i} = [pvalStr{i} ' ns'];
        end
        if p.Results.printP
            if round(tmpStats.("P-value"),4) == 0
                printedPval = ['p < 0.0001', pvalStr{i}];
            else
                printedPval = ['p = ',num2str(round(tmpStats.("P-value"),4)), pvalStr{i}];
            end
        else
            printedPval = pvalStr{i};
        end
        hold on
        plot([p.Results.comparison(i,1) p.Results.comparison(i,2)], repmat(distanceSignficancePlot(i),1,2),'k');
        text((p.Results.comparison(i,1)+p.Results.comparison(i,2))/2, (distanceSignficancePlot(i)),printedPval,...
            'HorizontalAlignment','Center','VerticalAlignment','top','FontName','Calibri');
    end
end

%---------------------------
% Finalizing plot
if p.Results.writeRatio == 1
    N = [];
    for i = 1:length(b)
        N = [N; b(i).YData];
    end
    N = sum(N,1);
    for i = 1:length(b)
        xtips = b(i).XEndPoints;
        ytips = b(i).YEndPoints;
        labels = string(round(ytips./N,3));
        text(xtips,ytips,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end
end
if ~isempty(p.Results.writeBarNames)
    for i = 1:length(p.Results.writeBarNames)
        xtips = b(i).XEndPoints;
        ytips = b(i).YEndPoints;
        labels = string(p.Results.writeBarNames{i});
        text(xtips,ytips,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom','FontSize',8)
    end
end
if ~isempty(p.Results.legend)
    legend(p.Results.legend);
    if ~isempty(p.Results.legendBox)
        if p.Results.legendBox == 0
            legend('boxoff');
        end
    end
    legend('location',p.Results.legendLoc);
end
ylabel(p.Results.ylabel);
title(p.Results.title);
xlim(p.Results.xlim);
if ~isempty(p.Results.ylim); ylim(p.Results.ylim); end
if strcmp(p.Results.grid, 'on'); grid on; end
if strcmp(p.Results.grid, 'off'); grid off; end
end

function [pVal, type, chi2stat, df] = computeChi2(X,alpha)
    % H0 : the variables are independent, there is no relationship between the two categorical variables. 
    % Knowing the value of one variable does not help to predict the value of the other variable
    %
    % H1 : the variables are dependent, there is a relationship between the two categorical variables. 
    % Knowing the value of one variable helps to predict the value of the other variable
    %
    % If the difference between the observed frequencies and the expected frequencies is small,
    % we cannot reject the null hypothesis of independence and thus we cannot reject the fact 
    % that the two variables are not related. On the other hand, if the difference between the 
    % observed frequencies and the expected frequencies is large, we can reject the null hypothesis
    % of independence and thus we can conclude that the two variables are related.
    
    if nargin == 2; alpha = 0.05; end
    if size(X,1) <= 1 || size(X,2) <= 1; error('Contingency table is missing values.'); end
    
    d = X; % Creating the table
    df = (size(d,1)-1)*(size(d,2)-1); % Computing the degrees of freedom; (n-row - 1) * (n-col - 1)
    
    N = sum(sum(d)); % Grand total
    ef = []; % Expected frequencies
    for i1 = 1:size(d,1)
        for i2 = 1:size(d,2)
            ef(i1,i2) = sum(d(i1,:),2)*sum(d(:,i2))/N;
        end
    end
    if sum(sum(ef < 5)) >= 1
        warning('Less than 5 samples in the expected frequency for one of the conditions. Switching to Fisher''s exact test.');
        type = 'Fishers exact test';
        [~,pVal,chi2stat] = fishertest(d);
    else
        type = 'Chi-squared test';
        chi2stat = []; % Chi-square table for each variable
        for i1 = 1:size(ef,1)
            for i2 = 1:size(ef,2)
                if size(ef,2) == 2; chi2stat(i1,i2) = ((abs(d(i1,i2)-ef(i1,i2))-0.5)^2)/ef(i1,i2); disp('Yates correction applied.'); end
                if size(ef,2) > 2; chi2stat(i1,i2) = ((d(i1,i2)-ef(i1,i2))^2)/ef(i1,i2); end
            end
        end
    
        chi2 = sum(sum(chi2stat)); % Pooled chi-square
        pVal = 1 - chi2cdf(chi2,df); % Computing the p-value
        if pVal == 0; pVal = chi2cdf(chi2,df,'upper'); end
    end
    
    
    if pVal < alpha; disp('Reject the null-hypothesis. The variables are dependent, there is a relationship between the two categorical variables.'); end
    if pVal >= alpha; disp('Cannot reject the null-hypothesis. The variables are independent, there is no relationship between the two categorical variables.'); end
end