function T_out = bbar_meanTable(T, groupVars)
% bbar_meanTable  Group table rows and take means of numeric columns.
%
%   T_out = bbar_meanTable(T, groupVars)
%
%   - T is a table.
%   - groupVars is a string scalar, string array, char, or cellstr of
%     variable names that define the groups (e.g. "ID", ["ID","Age"]).
%
%   For each unique combination of groupVars, this function:
%     * Computes the mean (omitnan) of all numeric-like columns.
%     * For string-like / categorical columns that are NOT in groupVars,
%       it requires that values are unique within each group.
%       If not unique, it throws an error.
%
%   Output:
%     - T_out is a table with one row per unique group, containing:
%         groupVars (as is)
%         all other variables aggregated appropriately.
%
%   Example:
%       % One row per ID, Age, with mean HeightDifference etc.
%       T_id = bbar_meanTable(T, ["ID", "Age"]);
%

    if ~istable(T)
        error('Input T must be a table.');
    end

    if nargin < 2 || isempty(groupVars)
        error('You must provide groupVars (grouping variable names).');
    end

    % Normalize groupVars to cellstr
    if isstring(groupVars)
        groupVars = cellstr(groupVars);
    elseif ischar(groupVars)
        groupVars = {groupVars};
    elseif iscellstr(groupVars)
        % ok
    else
        error('groupVars must be char, string, or cellstr of variable names.');
    end

    % Check that all grouping variables exist
    if ~all(ismember(groupVars, T.Properties.VariableNames))
        missing = setdiff(groupVars, T.Properties.VariableNames);
        error('These grouping variables are not in the table: %s', strjoin(missing, ', '));
    end

    % Get group indices and the table of unique group combinations
    [G, groupTable] = findgroups(T(:, groupVars));
    nGroups          = height(groupTable);

    % Start output with the grouping variables
    T_out = groupTable;

    allVars  = T.Properties.VariableNames;
    nonGroup = setdiff(allVars, groupVars, 'stable');

    % Loop over all non-grouping variables
    for i = 1:numel(nonGroup)
        vName = nonGroup{i};
        col   = T.(vName);

        % Prepare storage for aggregated column
        aggCol = [];

        % ---- Numeric-like variables: mean per group ----
        if isnumeric(col) || islogical(col) || isduration(col) || isdatetime(col)

            aggCol = splitapply(@(x) mean(x, 'omitnan'), col, G);

        % ---- String-like variables: must be unique per group ----
        elseif isstring(col)
            aggCol = strings(nGroups, 1);
            for g = 1:nGroups
                idx = (G == g);
                vals = col(idx);
                u = unique(vals);
                if numel(u) ~= 1
                    error('Variable "%s" has multiple string values within group %d.', vName, g);
                end
                aggCol(g) = u(1);
            end

        elseif iscellstr(col)
            aggCol = cell(nGroups, 1);
            for g = 1:nGroups
                idx = (G == g);
                vals = col(idx);
                u = unique(vals);
                if numel(u) ~= 1
                    error('Variable "%s" has multiple cellstr values within group %d.', vName, g);
                end
                aggCol{g} = u{1};
            end

        elseif iscategorical(col)
            aggCol = categorical(zeros(nGroups,1)); % temp init
            for g = 1:nGroups
                idx = (G == g);
                vals = col(idx);
                u = unique(vals);
                if numel(u) ~= 1
                    error('Variable "%s" has multiple categorical values within group %d.', vName, g);
                end
                aggCol(g) = u(1);
            end

        elseif ischar(col)
            % char array column (rare; treat as row-wise strings)
            % Convert to cellstr for convenience
            colCell = cellstr(col);
            aggCol  = cell(nGroups, 1);
            for g = 1:nGroups
                idx = (G == g);
                vals = colCell(idx);
                u = unique(vals);
                if numel(u) ~= 1
                    error('Variable "%s" has multiple char values within group %d.', vName, g);
                end
                aggCol{g} = u{1};
            end
            aggCol = string(aggCol); % store as string

        else
            error('Variable "%s" has unsupported type (%s).', vName, class(col));
        end

        % Attach aggregated column to output table
        T_out.(vName) = aggCol;
    end
end