% bbar_raincloud()    - Creates a raincloud plot.
%                         
%
% Usage:
%   >> [cmap, stats] = bbar_raincloud(data, varargin)
%   >> data{1} = ran
%   >> [cmap, stats] = bbar_raincloud(data, 'showMean', 1, 'writeMean', 1, 'stats', 'ttest')
%
% Inputs:
%   data                    - (Required) Data in cells:
%                             data{1} = 1x30 cells <- Condition 1
%                             data{2} = 1x30 cells <- Condition 2
%   writeMean               - Plots text that displays the mean (x-bar) (default: 0)
%   writeMedian             - Plots text that displays the median (x-tilde) (default: 0)
%   BoxPlot                 - Adds a box plot (default: 1)
%   BoxFaceAlpha            - Box plot alpha level (default: .5)
%   BoxWidth                - Width of boxes in box plot (default: .25)
%   Scatter                 - Adds a scatter plot (default: 1)
%   Normal                  - Adds a ksdensity plot (default: 1)
%   BarPlot                - Adds a barplot (default: 0)
%   Alpha                   - Alpha level of plots (default: .5)
%   jitterAmount            - Determines jitter level of scatterplot (default: .1)
%   cmap                    - Determines the colormap of the plots (default: bbar_makecmap())
%   showMean                - Plots a line through the mean of conditions (default: 0)
%   showMedian              - Plots a line through the median of conditions (default: 0)
%   showInteraction         - Plots a line through through conditions of each participant (default: 0)
%   nexttile                - In case of using 'tiledlayout', set this to
%                             true/1 so that it plots it as 'nexttile'
%   Bandwidth               - Bandwidth of the ksdensity plot (default: [])
%   xticks                  - Set the ticks on xline (default: based on data length)
%   xtickslabels            - Set the labels, eg. {'', 'Pre', 'Post'}
%   xlim                    - Set the extremes of the x axis.
%   plotLsline              - Plots the least-square line to show a mean
%                             (default: 0)
%   stats                   - Compute the statistics, either 'wilcoxon',
%                             'steeldwass', 'ttest' of 'chi2' (default: [])
%                             if chi2 is used, make sure to supply two
%                             lines data input, e.g. 2 x 1 cell array:
%                             input = [{8},{5};{4},{7}]. Assuming here that the
%                             first row are the Incorrect responses.
%   comparison              - Comparison matrix, eg. [1 2; 2 3; 1 3] which compares
%                             condition 1 with 2, and 2 with 3, and 1 with 3 (default: [])
%   multcmp                 - If there's more than one comparison, it will
%                             automatically correct for multiple
%                             comparison. Determine the method:
%                             'holm', 'hochberg', 'hommel', 
%                             'bonferroni', 'BH', 'BY', 'fdr'(default), 'sidak' or 'none'.
%                             Number of comparisons can be custom if used
%                             like a cell: {'BY',3}: Use 'BY', correct for
%                             3 comparisons.
%   rmOutlier               - Remove outliers. Set the z-score difference.
%                             (default: 0). If non-zero, outliers are
%                             typically set to 3.
%   paired                  - Determine if 'paired' (default: 0)
%   printP                  - Determine if exact p-values should be printed
%                             with significance stars (default: 0)
%   ylabel                  - The unit of data (default: 'magnitude')
%   xlabel                  - The unit-change of data (default: '')
%   title                   - Add title to plot (default: '')
%   fontName                - Set the font of the plots (default: 'Arial')
%   grid                    - Add grid to plot (default: 'on')
%   bgColor                 - Set background color of figure (default: 'w')
%   StatsYPos               - The location of the statistics in the plot.
%                             Must include two values. First one determines
%                             the y-position. Second one determines the 
%                             distance between the conditions (if multiple), eg. [3 0.1].
%   fPosition               - Determine the position of the figure. This
%                             affects where the multcomp plots,
%                             e.g., ("p = 0.0001"), will appear.
%
%
% Outputs:
%   cmap                    - The color map used in the plot
%   tmpStats                - The statistical details, if computed
%
%
% Author: Zakaria Djebbara, 2022
%
%

function [cmap, tmpStats] = bbar_raincloud(data, varargin)

defaultXTicksLabels{1} = '';
for i = 0:length(data)
    defaultXTicksLabels{1,i+1} = ['Condition ' num2str(i)];
end
defaultXTicksLabels{end+1} = '';
defaultXTicks = 0:length(data)+1;


p = inputParser;
addRequired(p, 'data', @iscell);
addOptional(p, 'rmOutlier', 0, @isnumeric);
addOptional(p, 'writeMean', 0, @isnumeric);
addOptional(p, 'writeMedian', 0, @isnumeric);
addOptional(p, 'BoxPlot', 1, @isnumeric);
addOptional(p, 'BoxFaceAlpha', .5, @isnumeric);
addOptional(p, 'BoxWidth', .25, @isnumeric);
addOptional(p, 'Scatter', 1, @isnumeric);
addOptional(p, 'ScatterSize', 15, @isnumeric);
addOptional(p, 'Normal', 1, @isnumeric);
addOptional(p, 'BarPlot', 0, @isnumeric);
addOptional(p, 'alpha', .5, @isnumeric);
addOptional(p, 'jitterAmount', .1, @isnumeric);
addOptional(p, 'cmap', bbar_makecmap([51/255, 92/255, 103/255; 224/255, 159/255, 62/255; 84/255, 11/255, 14/255], size(data,2)),...
    @(x) (size(data,2) <= size(x,1)) && isnumeric(x));
addOptional(p, 'showMean', 0);
addOptional(p, 'showMedian', 0);
addOptional(p, 'showInteraction', 0);
addOptional(p, 'nexttile', 0);
addOptional(p, 'Bandwidth', [], @isnumeric);
addOptional(p, 'xticks', defaultXTicks, @isnumeric);
addOptional(p, 'xtickslabels', defaultXTicksLabels, @iscell);
addOptional(p, 'writeXTicks', 1, @isnumeric);
addOptional(p, 'xlim', [.5 length(data)+.75], @isnumeric);
addOptional(p, 'ylim', [], @isnumeric);
addOptional(p, 'plotLsline', 0, @isnumeric);
addOptional(p, 'stats', '', @ischar);
addOptional(p, 'comparison', [], @isnumeric);
addOptional(p, 'multcmp','fdr', @(x) ischar(x) || iscell(x));
addOptional(p, 'paired', 0, @isnumeric);
addOptional(p, 'ylabel', 'magnitude', @ischar);
addOptional(p, 'xlabel', '', @ischar);
addOptional(p, 'title', '', @ischar);
addOptional(p, 'fontName', 'Arial', @ischar);
addOptional(p, 'grid', 'on');
addOptional(p, 'FDR', 0.05, @isnumeric);
addOptional(p, 'bgColor', 'w');
addOptional(p, 'printP',0,@isnumeric);
addOptional(p, 'StatsYPos',[],@(x) (sum(size(x)) == 3) && isnumeric(x));
addOptional(p, 'fPosition',[200 200 600 400],@isvector);
parse(p,data,varargin{:});

% Creating a default colormap
cmap = p.Results.cmap;

% Making sure the 'nexttile' is an option
if strcmp(p.Results.nexttile, 'on')
    nexttile;
elseif p.Results.nexttile ~= 0
    nexttile(p.Results.nexttile);
elseif p.Results.nexttile == 0
    figure('Position',p.Results.fPosition)
end

% Making sure the structure is "tall"/columns
for i = 1:length(data)
    % In case of quadratic sizes, give warning
    if size(data{i},1) == size(data{i},2)
        warning('Make sure the data is structured so that each column is a condition, and each row a measure.');
    end
    if size(data{i},1) < size(data{i},2)
        data{i} = data{i}';
    end
end

% If remove outlier is on
if p.Results.rmOutlier
    for i = 1:length(data)
        data{i} = rmoutliers(data{i},"quartiles");
    end
end 

% If not equally long, fill in with NaN's
tmpList = [];
for i = 1:length(data)
    tmpList = [tmpList; length(data{i})];
end
if min(tmpList) ~= max(tmpList)
    for i = 1:length(data)
        data{i}(length(data{i}):max(tmpList)) = NaN;
    end
end

% Making sure the cells are dispatched correctly
if ~strcmp(p.Results.stats,'chi2')
    try 
        cell2mat(data{i})
    catch
        for i = 1:length(data)
            data{i} = mat2cell(data{i},repmat(1,1,length(data{i})));
        end
    end
else
    % if size(data{i},2) < 2
    %     error('When using Chi2, please make sure the data is structured so that data{1} = [25, 50, 90; 75, 50, 10] where the first row represent the Hits, and second row represent the Miss. Each column is a participant.');
    % end
    
    for i2 = 1:length(data)
        data{i2} = [mat2cell(data{i2}(:,1),repmat(1,1,length(data{i2}))), mat2cell(data{i2}(:,2),repmat(1,1,length(data{i2})))];
    end

end

% Generating the jitter-pattern
jitter = remapnumbers(rand(1,length(cell2mat(data{i}))), -p.Results.jitterAmount, p.Results.jitterAmount);

% Making sure that if 'stats' has been declared, the comparison should too.
if ~isempty(p.Results.stats) && isempty(p.Results.comparison); error('Please provide a comparison if you wish to compute the stats. E.g. [1 2; 2 3]'); end
if isempty(p.Results.stats) && ~isempty(p.Results.comparison); error('Please provide a statistical test if you wish to compute the stats. E.g. wilcoxon.'); end

% Making sure the multcmp is output correctly
if iscell(p.Results.multcmp)
    for i = 1:length(p.Results.multcmp)
        if ~ischar(p.Results.multcmp{1}) || ~isnumeric(p.Results.multcmp{2})
            error('Error in multcmp: first provide the method, and then the number of DF.');
        end
        
        DFnum = p.Results.multcmp{2};
        multcmp = p.Results.multcmp{1};

    end
else
    multcmp = p.Results.multcmp;
end

% Box Plot
if p.Results.BoxPlot 
    tmpBoxData = [];
    for i = 1:length(data)
        tmpBoxData = [tmpBoxData, cell2mat(data{i})];
    end
    g = repmat({''}, size(tmpBoxData,1)*size(tmpBoxData,2), 1);
    xBox = (1:size(tmpBoxData,2)) - .1;
    for i = 1:size(tmpBoxData,2)
        b = boxchart( repmat(xBox(i),size(tmpBoxData,1),1), tmpBoxData(:,i), 'MarkerStyle','none',...
            'BoxFaceColor',p.Results.cmap(i,:), 'BoxWidth', p.Results.BoxWidth, 'BoxFaceAlpha', p.Results.alpha);
        hold on;
    end
    set(gca,'xticklabel',g)
end

% Interaction plot
% Assuming each row represents a participant
if p.Results.showInteraction
    placeholder = [];
    xPlaceholder = [];
    for i = 1:length(data)
        placeholder(:,i) = cell2mat(data{i}(:));
        xPlaceholder = [xPlaceholder; i + jitter];
    end
    for i = 1:length(cell2mat(data{1}(:)))
        if p.Results.Scatter
            plot(xPlaceholder(:,i)',[placeholder(i,:)],'LineWidth',0.15, 'Color', [0.5 0.5 0.5]);
        else
            plot([1:length(data)],[placeholder(i,:)],'LineWidth',0.15, 'Color', [0.5 0.5 0.5]);
        end
    end
end

% Scatter Plot
if p.Results.Scatter
    for i = 1:length(data)
        xScatter{i} = [i+jitter(1:length(cell2mat(data{i})));];
        if strcmp(p.Results.stats,'chi2')
            scatter(xScatter{i},...
                [cell2mat(data{i}(:,1))],...
                p.Results.ScatterSize,p.Results.cmap(i,:),'filled','MarkerFaceAlpha',p.Results.alpha);
        else 
            scatter(xScatter{i},...
                [cell2mat(data{i})],...
                p.Results.ScatterSize,p.Results.cmap(i,:),'filled','MarkerFaceAlpha',p.Results.alpha);
        end
        hold on;
    end
end

% linear fit (it's working fine, but needs to be reworked so that one can
% choose between fitting within one condition, or across all conditions
if p.Results.plotLsline
    try 
        h = lsline;
    catch
        error('Using the linear regression requires setting ''scatter'' to 1');
    end
    for i1 = 1:length(h)
        % h(i1).Color = mycmap(i1,:);
        h(i1).Color = 'k';
        h(i1).LineStyle = ':';
        h(i1).LineWidth = 1;
    end
end

% Bar Plot
if p.Results.BarPlot
    tempData = [];
    tempCmap = [];
    for i = 1:length(data)
        tempData = [tempData, mean(cell2mat(data{i}(:,1)),'omitnan')];
        tempCmap = [tempCmap; p.Results.cmap(i,:)];     
    end
    bb = bar(tempData, 'barwidth',p.Results.BoxWidth,'FaceColor','flat','EdgeColor','none', 'FaceAlpha', p.Results.alpha);
    bb.CData = tempCmap;
    bb.BaseLine.LineStyle = ":";
    hold on;
end

% Density (bell-curve) plot
if p.Results.Normal
    for i = 1:length(data)
        [f, xi] = ksdensity(cell2mat(data{i}),'Bandwidth', p.Results.Bandwidth);
        pfig = patch(remapnumbers(f, 0, .5) + i, xi, p.Results.cmap(i,:), 'FaceAlpha',p.Results.alpha);
    end
end

% Computing the mean for each condition
ax1 = gca;
tempmean = [];
tempmedian = [];
tempmin = [];
if strcmp(p.Results.stats, 'chi2')
    for i = 1:length(data)
        tempmean = [tempmean; mean(cell2mat(data{i}(1,:)),'omitnan')];
        tempmedian = [tempmedian; median(cell2mat(data{i}(1,:)),'omitnan')];
        tempmin = [tempmin; min(cell2mat(data{i}(1,:)))];
    end
else
    for i = 1:length(data)
        tempmean = [tempmean; mean(cell2mat(data{i}),'omitnan')];
        tempmedian = [tempmedian; median(cell2mat(data{i}),'omitnan')];
        tempmin = [tempmin; min(cell2mat(data{i}))];
    end
end

% Display mean values
if p.Results.showMean
    if ~isempty(p.Results.comparison)
        temp_comparison = cell(size(p.Results.comparison,1),1);
        if size(temp_comparison,1) == 1 && size(temp_comparison,2) == 1
            plot((1:length(tempmean))-(p.Results.BoxWidth/2),tempmean,'k','LineWidth', 1.5);
        else
            for i = 1:size(p.Results.comparison,1)
                temp_comparison{i} = p.Results.comparison(i,1):p.Results.comparison(i,2);
            end
            adj_C = test_similarity(temp_comparison);
            adj_C = adj_C(~cellfun('isempty',adj_C));
            for i = 1:size(adj_C,1)
                plot(adj_C{i}-(p.Results.BoxWidth/2),tempmean(adj_C{i}),'k','LineWidth', 1.5);
            end
        end
    else
        plot((1:length(tempmean))-(p.Results.BoxWidth/2),tempmean,'k','LineWidth', 1.5);
    end
end
% Display median values
if p.Results.showMedian
    if ~isempty(p.Results.comparison)
        temp_comparison = cell(size(p.Results.comparison,1),1);
        if size(temp_comparison,1) == 1 && size(temp_comparison,2) == 1
            plot((1:length(tempmedian))-(p.Results.BoxWidth/2),tempmedian,':k','LineWidth', 1.5);
        else
            for i = 1:size(p.Results.comparison,1)
                temp_comparison{i} = p.Results.comparison(i,1):p.Results.comparison(i,2);
            end
            adj_C = test_similarity(temp_comparison);
            adj_C = adj_C(~cellfun('isempty',adj_C));
            for i = 1:size(adj_C,1)
                plot(adj_C{i}-(p.Results.BoxWidth/2),tempmedian(adj_C{i}),':k','LineWidth', 1.5);
            end
        end
    else
        plot((1:length(tempmedian))-(p.Results.BoxWidth/2),tempmedian,':k','LineWidth', 1.5);
    end
end
% Write the mean values
if p.Results.writeMean
    for i = 1:length(tempmean)
        q = quantile(cell2mat(data{i}(1,:)),3);
        I = find(cell2mat(data{i}(1,:)) <= q(1)-1.5*(q(3)-q(1)));
        if isempty(I)
            I = min(cell2mat(data{i}(1,:)));
        elseif ~isempty(I)
            tmp = cell2mat(data{i}(1,:));
            tmp(I) = [];
            I = min(tmp);
        end
        text(i-(p.Results.BoxWidth/2), ax1.YTick(end),...
            ['$\bar{x}$:' num2str(round(tempmean(i),2))],'interpreter','latex',...
            'VerticalAlignment','top','HorizontalAlignment','center','FontSize',12);
%         ht = text(i-(p.Results.BoxWidth/2), tempmedian(i),...
%             ['$\bar{x}$: ' num2str(round(tempmean(i),2))],'interpreter','latex',...
%             'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',14);
%         set(ht,'Rotation',90);
    end
end
% Write the median values
if p.Results.writeMedian
    for i = 1:length(tempmean)
        q = quantile(cell2mat(data{i}),3);
        I = find(cell2mat(data{i}) <= q(1)-1.5*(q(3)-q(1)));
        if isempty(I)
            I = min(cell2mat(data{i}));
        elseif ~isempty(I)
            tmp = cell2mat(data{i});
            tmp(I) = [];
            I = min(tmp);
        end
        text(i-(p.Results.BoxWidth/2), I,...
            ['$\tilde{x}$: ' num2str(round(tempmean(i),2))],'interpreter','latex',...
            'VerticalAlignment','top','HorizontalAlignment','center','FontSize',12);
    end
end

if isempty(p.Results.StatsYPos)
    distanceSignficancePlot = [];
    for i = 1:size(p.Results.comparison,1)
        distanceSignficancePlot = [distanceSignficancePlot; ax1.YTick(end)+(ax1.YTick(end)-ax1.YTick(end-1))*i];
    end
else
    distanceSignficancePlot = [];
    for i = 1:size(p.Results.comparison,1)
        if i == 1; distanceSignficancePlot = [distanceSignficancePlot; p.Results.StatsYPos(1)]; end
        if i ~= 1; distanceSignficancePlot = [distanceSignficancePlot; distanceSignficancePlot(end)+p.Results.StatsYPos(2)]; end
    end
end

% Stats
% --------------------------------------------------
%
% TO DO:
% Linear regression 
%
tmpStats{7} = [];
pvalStr = cell(size(p.Results.comparison,1),1);
condNames = p.Results.xtickslabels;

% Chi2
% ------------------
if strcmp(p.Results.stats, 'chi2')
    for i = 1:size(p.Results.comparison,1)
        % Check if input is binary
        tmpD1 = cell2mat(data{p.Results.comparison(i,1)});
        tmpD2 = cell2mat(data{p.Results.comparison(i,2)});
        tmpD1 = ~isnan(tmpD1);
        tmpD2 = ~isnan(tmpD2);
        if size(data{1},1) ~= 2
            error('Please make sure to input counts as cells, e.g., 2 x 1 cell array: input = [{8};{4}]');
        end
        % Computing Chi2
        [pVal, type, chi2stat] = computeChi2(cell2mat([data{p.Results.comparison(i,1)},data{p.Results.comparison(i,2)}]),p.Results.alpha);
        pVal = round(pVal,4,'decimal');
        tmpStats{i,1} = type;
        tmpStats{i,2} = [tmpStats{i,2}; condNames(p.Results.comparison(i,1)) + " vs " + condNames(p.Results.comparison(i,2))];
        if strcmp(type, 'Chi-squared test'); tmpStats{i,3} = [tmpStats{i,3}; sum(sum(chi2stat))]; tci = 'NA'; end
        if strcmp(type, 'Fishers exact test'); tmpStats{i,3} = [tmpStats{i,3}; chi2stat.OddsRatio]; tci = chi2stat.ConfidenceInterval; end
        tmpStats{i,4} = [tmpStats{i,4}; [round(tempmean(p.Results.comparison(i,1)),2) round(tempmean(p.Results.comparison(i,2)),2)]];
        tmpStats{i,5} = [tmpStats{i,5}; pVal];
    end
    
    tmpStats = cell2table(tmpStats);
    tmpStats = renamevars(tmpStats,["tmpStats1","tmpStats2","tmpStats3","tmpStats4","tmpStats5","tmpStats6"], ...
    ["Stats test","Comparison","Chi/Odds ratio","Means","P-value","Cohens d"]);
end

% Steel-Dwass 
% ------------------
% Wilcoxon compares medians, this one compares means.
if strcmp(p.Results.stats,'steel') || strcmp(p.Results.stats,'dwass') || strcmp(p.Results.stats,'steeldwass')
    disp('Working on this one still. Functions are done but implementation is missing.')

    
        % if p.Results.paired == 0
        %     preparedObservations = [];
        %     preparedGroups = [];
        %     for i = 1:size(p.Results.xtickslabels,2)
        %         preparedGroups = [preparedGroups; repmat(string(p.Results.xtickslabels{i}),length(cell2mat(data{i})),1)];
        %         preparedObservations = [preparedObservations; cell2mat(data{i})];
        %     end
        %     assert(length(preparedGroups) == length(preparedObservations));
        %     [pVal, stats] = steel_dwass_unpaired(preparedObservations, preparedGroups);
        %     tmpStats{i,1} = ["Steel-Wass Unpaired test"];
        % elseif p.Results.paired == 1
        %     preparedObservations = [];
        %     preparedGroups = [];
        %     for i = 1:size(p.Results.xtickslabels,2)
        %         preparedGroups = [preparedGroups; repmat(string(p.Results.xtickslabels{i}),length(cell2mat(data{i})),1)];
        %         preparedObservations = [preparedObservations; cell2mat(data{i})];
        %     end
        %     assert(length(preparedGroups) == length(preparedObservations));
        %     [pVal, stats] = steel_dwass_paired(preparedObservations, preparedGroups);
        %     tmpStats{i,1} = ["Steel-Wass Unpaired test"];
        % end
    
end

% Wilcoxon
% ------------------
if strcmp(p.Results.stats, 'wilcoxon')
    for i = 1:size(p.Results.comparison,1)
        % Wilcoxon
        if p.Results.paired == 0
            [pVal,h,stats] = ranksum(cell2mat(data{p.Results.comparison(i,1)}), cell2mat(data{p.Results.comparison(i,2)}),'method','approximate','tail','both');
            tmpStats{i,1} = ["Wilcoxon Unpaired test"];
        else
            [pVal,h,stats] = signrank(cell2mat(data{p.Results.comparison(i,1)}), cell2mat(data{p.Results.comparison(i,2)}),'method','approximate','tail','both');
            tmpStats{i,1} = ["Wilcoxon Paired test"];
        end
        pVal = round(pVal,4,'decimal');
        tmpStats{i,2} = [tmpStats{i,2}; condNames(p.Results.comparison(i,1)) + " vs " + condNames(p.Results.comparison(i,2))];
        if ~isfield(stats,'zval'); stats.zval = NaN; end
        tmpStats{i,3} = [tmpStats{i,3}; stats.zval];
        tmpStats{i,4} = [tmpStats{i,4}; [round(tempmean(p.Results.comparison(i,1)),2) round(tempmean(p.Results.comparison(i,2)),2)]];
        % tmpPval{1} = [tmpPval{1}; condNames(p.Results.comparison(i,1)) + " vs " + condNames(p.Results.comparison(i,2))];
        tmpStats{i,5} = [tmpStats{i,5}; pVal];
        [d, ci_lower, ci_upper] = computeCohensD(cell2mat(data{p.Results.comparison(i,1)}), cell2mat(data{p.Results.comparison(i,2)}));
        tmpStats{i,6} = [tmpStats{i,6}; d];
        tmpStats{i,7} = [tmpStats{i,7}; [ci_lower, ci_upper]];
    end
    tci = 'NA';
    tmpStats = cell2table(tmpStats);
    tmpStats = renamevars(tmpStats,["tmpStats1","tmpStats2","tmpStats3","tmpStats4","tmpStats5","tmpStats6","tmpStats7"], ...
        ["Stats test","Comparison","Z-Val","Means","P-value","Cohens d","CI95%"]);
end

% T-Test 
% ------------------
if strcmp(p.Results.stats, 'ttest')
    tmpStats{9} = [];
    for i = 1:size(p.Results.comparison,1)
        % Test for variance (Levenes)
        LevP = vartestn([cell2mat(data{p.Results.comparison(i,1)}), cell2mat(data{p.Results.comparison(i,2)})],...
                'TestType','LeveneAbsolute','Display','off');
        if LevP <= 0.05; warning(['Comparison [' num2str(p.Results.comparison(i,1)) ' ' num2str(p.Results.comparison(i,2)) ...
                '] did not pass Levenes variance test (p = ' num2str(LevP) '. Data has different variance. Consider using ''Wilcoxon'' test --> ''stats'',''wilcoxon''.']); end
        % Test for normality
        [~, pLillie] = lillietest(cell2mat(data{p.Results.comparison(i,1)}));
        if pLillie > 0.05; warning(['Data [' num2str(p.Results.comparison(i,:)) '] is normal (Lilliefors p = ' num2str(pLillie) ').']); end
        if pLillie <= 0.05; warning(['Data [' num2str(p.Results.comparison(i,:)) '] did not pass Lilliefors test.'...
                ' Data is not normal (Lilliefors p = ' num2str(pLillie) '). Consider using ''Wilcoxon'' test --> ''stats'',''wilcoxon''.']); end
        % T-Test
        if p.Results.paired == 0
            [h, pVal, tci, stats] = ttest2(cell2mat(data{p.Results.comparison(i,1)}), cell2mat(data{p.Results.comparison(i,2)}));
            tmpStats{i,1} = 'Unpaired t-test';
        else 
            [h, pVal, tci, stats] = ttest(cell2mat(data{p.Results.comparison(i,1)}), cell2mat(data{p.Results.comparison(i,2)}));
            tmpStats{i,1} = 'Paired t-test';
        end
            pVal = round(pVal,4,'decimal');
            tmpStats{i,2} = [tmpStats{i,2}; condNames(p.Results.comparison(i,1)) + " vs " + condNames(p.Results.comparison(i,2))];
            tmpStats{i,3} = [tmpStats{i,3}; stats.sd];
            tmpStats{i,4} = [tmpStats{i,4}; tci'];
            % tmpPval{1} = [tmpPval{1}; condNames(p.Results.comparison(i,1)) + " vs " + condNames(p.Results.comparison(i,2))];
            tmpStats{i,5} = [tmpStats{i,5}; pVal];
            tmpStats{i,6} = [tmpStats{i,6}; LevP];
            tmpStats{i,7} = [tmpStats{i,7}; pLillie];
            [d, ci_lower, ci_upper] = computeCohensD(cell2mat(data{p.Results.comparison(i,1)}), cell2mat(data{p.Results.comparison(i,2)}));
            tmpStats{i,8} = [tmpStats{i,8}; d];
            tmpStats{i,9} = [tmpStats{i,9}; [ci_lower, ci_upper]];
    end
    tmpStats = cell2table(tmpStats);
    tmpStats = renamevars(tmpStats,["tmpStats1","tmpStats2","tmpStats3","tmpStats4","tmpStats5","tmpStats6","tmpStats7","tmpStats8","tmpStats9"], ...
        ["Stats test","Comparison","SD","TCI","P-value","Levene's P-value","Lilliefors' P-value", "Cohen's d","CI95%"]);
end

% Levene Variance Test
% -------------------
if strcmp(p.Results.stats, 'levene')
    tmpStats{4} = [];
    for i = 1:size(p.Results.comparison,1)
        [pVal,stats] = vartestn([cell2mat(data{p.Results.comparison(i,1)}), cell2mat(data{p.Results.comparison(i,2)})],...
                'TestType','LeveneAbsolute','Display','off');
        pVal = round(pVal,4,'decimal');
        warning('Computing and displaying Levene Variance Test.');
    end
    tmpStats{2} = [tmpStats{2}; condNames(p.Results.comparison(i,1)) + " vs " + condNames(p.Results.comparison(i,2))];
    tmpStats{3} = [tmpStats{3}; stats.fstat];
    tmpStats{4} = [tmpStats{4}; stats.df];
    tmpPval{1} = [tmpPval{1}; condNames(p.Results.comparison(i,1)) + " vs " + condNames(p.Results.comparison(i,2))];
    tmpPval{2} = [tmpPval{2}; pVal];
end

% Correcting for multiple comparison if comparions are more than one.
% -------------------
if ~isempty(p.Results.stats)
    if strcmp(multcmp,'none')
        adj_p = tmpStats.("P-value");
    else
        if exist('DFnum','var')
            adj_p = pval_adjust(tmpStats.("P-value"),multcmp, DFnum);
        else
            adj_p = pval_adjust(tmpStats.("P-value"),multcmp);
        end
        if size(p.Results.comparison,1) > 1
            warning(['p-values are corrected for multiple comparison using ' multcmp '.']);
        end
    end
    tmpStats.("P-value") = adj_p;
    for i = 1:size(p.Results.comparison,1)
        if tmpStats.("P-value")(i) <= 0.1 && tmpStats.("P-value")(i) > 0.05
            pvalStr{i} = ['' pvalStr{i} ' ·'];
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
                printedPval = ['p = ', num2str(round(tmpStats.("P-value")(i),4),'%.4f'), pvalStr{i}];
            end
        else
            printedPval = pvalStr{i};
        end
        plot([p.Results.comparison(i,1) p.Results.comparison(i,2)], repmat(distanceSignficancePlot(i),1,2),'k');
        text((p.Results.comparison(i,1)+p.Results.comparison(i,2))/2, (distanceSignficancePlot(i)),printedPval,...
            'HorizontalAlignment','Center','VerticalAlignment','top','FontName','Calibri','FontSize',8);
    end
end

% TODO
% array2table(tmpStats,"VariableNames", ...
%    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

% Finalizing plot
if p.Results.writeXTicks
    xticklabels([p.Results.xtickslabels]);
    xticks(p.Results.xticks+1);
end
ylabel(p.Results.ylabel);
xlabel(p.Results.xlabel);
title(p.Results.title);
if ~isempty(p.Results.xlim); xlim(p.Results.xlim); end
if ~isempty(p.Results.ylim); ylim(p.Results.ylim); end
% if p.Results.writeMean
%     y1 = ylim;
% %     ylim([min(pfig.Vertices(:,2))*4 y1(2)]);
%     ylim([min(tempmin)-(std(tempmin)*3) y1(2)]);
% end
if ~isempty(p.Results.ylim); ylim(p.Results.ylim); end
if strcmp(p.Results.grid, 'on'); grid on; end
if strcmp(p.Results.grid, 'off'); grid off; end

set(gcf,'Color',p.Results.bgColor);
set(gca,'fontname',p.Results.fontName);

% ax1.YTick = ax1.YTick(end)+(ax1.YTick(end)-ax1.YTick(end-1));
end

function adj_C = test_similarity(C)
    t = perms(1:size(C,1));
    t = t(:,1:2);
    for i1 = 1:size(t,1)
        [lia, lib] = ismember(C{t(i1,1)},C{t(i1,2)});
        if sum(lia) >= 2
            C{t(i1,1)} = unique([C{t(i1,1)}, C{t(i1,2)}]);
            C{t(i1,2)} = [];
        end
    end
    adj_C = C;
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

function [d, ci_lower, ci_upper] = computeCohensD(x1, x2, alpha)
    % Compute Cohen's d effect size with confidence intervals
    % x1: vector of observations for group 1
    % x2: vector of observations for group 2
    % alpha: significance level for confidence intervals (default is 0.05)
    
    if nargin < 3
        alpha = 0.05; % default significance level
    end
    
    % Compute means
    mean1 = mean(x1,'omitnan');
    mean2 = mean(x2,'omitnan');
    
    % Compute pooled standard deviation
    n1 = length(x1);
    n2 = length(x2);
    pooledStd = sqrt(((n1-1)*std(x1,'omitnan')^2 + (n2-1)*std(x2,'omitnan')^2) / (n1 + n2 - 2));
    
    % Compute Cohen's d
    d = (mean1 - mean2) / pooledStd;
    
    % Compute standard error of d
    sed = pooledStd * sqrt(1/n1 + 1/n2);
    
    % Compute confidence intervals
    ci_lower = d - norminv(1-alpha/2) * sed;
    ci_upper = d + norminv(1-alpha/2) * sed;
end

function [p_values, stats] = steel_dwass_unpaired(data, groups)
% STEEL_DWASS_UNPAIRED Performs Steel-Dwass test for unpaired data
% Returns uncorrected p-values and test statistics
% 
% Inputs:
%   data: Column vector of observations
%   groups: Column vector of group labels (same length as data)
%
% Outputs:
%   p_values: Vector of uncorrected p-values for each comparison
%   stats: Vector of test statistics for each comparison
%   comparison_labels: Cell array describing each comparison

unique_groups = unique(groups);
n_groups = length(unique_groups);
n_comparisons = nchoosek(n_groups, 2);
p_values = zeros(n_comparisons, 1);
stats = zeros(n_comparisons, 1);
comparison_labels = cell(n_comparisons, 1);

comparison = 1;
for i = 1:n_groups
   for j = (i+1):n_groups
       % Get data for groups i and j
       group1 = data(groups == unique_groups(i));
       group2 = data(groups == unique_groups(j));
       
       n1 = length(group1);
       n2 = length(group2);
       
       % Combine and rank
       combined = [group1; group2];
       [~, ranks] = sort(combined);
       rank_values = zeros(size(combined));
       rank_values(ranks) = 1:length(combined);
       
       % Calculate rank sum and test statistic
       R1 = sum(rank_values(1:n1));
       S = (R1 - n1*(n1+n2+1)/2) / sqrt(n1*n2*(n1+n2+1)/12);
       
       % Store test statistic and p-value
       stats(comparison) = S;
       p_values(comparison) = 2 * (1 - normcdf(abs(S)));
       comparison_labels{comparison} = [unique_groups(i) + " vs " + unique_groups(j)];
       
       comparison = comparison + 1;
   end
end
end

function [p_values, stats] = steel_dwass_paired(data, groups)
% STEEL_DWASS_PAIRED Performs Steel-Dwass test for paired data
% Returns uncorrected p-values and test statistics
% 
% Inputs:
%   data: Matrix where each column is a condition and each row is a subject
%   groups: Vector of group labels (length = number of columns in data)
%
% Outputs:
%   p_values: Vector of uncorrected p-values for each comparison
%   stats: Vector of test statistics for each comparison
%   comparison_labels: Cell array describing each comparison


unique_groups = unique(groups);
n_groups = length(unique_groups);
n_comparisons = nchoosek(n_groups, 2);
p_values = zeros(n_comparisons, 1);
stats = zeros(n_comparisons, 1);
comparison_labels = cell(n_comparisons, 1);

comparison = 1;
for i = 1:n_groups
   for j = (i+1):n_groups
       % Get paired differences
       diff = data(:,i) - data(:,j);
       
       % Remove zeros and get signs
       nonzero_diff = diff(diff ~= 0);
       signs = sign(nonzero_diff);
       
       % Rank absolute differences
       [~, ranks] = sort(abs(nonzero_diff));
       rank_values = zeros(size(nonzero_diff));
       rank_values(ranks) = 1:length(nonzero_diff);
       
       % Calculate signed rank sum
       W = sum(signs .* rank_values);
       
       % Calculate test statistic
       n = length(nonzero_diff);
       S = W / sqrt(n*(n+1)*(2*n+1)/6);
       
       % Store test statistic and p-value
       stats(comparison) = S;
       p_values(comparison) = 2 * (1 - normcdf(abs(S)));
       comparison_labels{comparison} = [unique_groups(i) + " vs " + unique_groups(j)];
       
       comparison = comparison + 1;
   end
end
end