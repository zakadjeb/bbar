function cell_out = bbar_meanCell(cell_in,direction)
% Default
if nargin == 1; direction = 1; end

% Check cell sizes
sizeMatrix = cellfun(@size, cell_in, 'UniformOutput', false);
for i = 1:length(sizeMatrix)-1
    if ~isequal(sizeMatrix{i}, sizeMatrix{i+1})
        error('All cells must have the same size.');
    end
end

% Stack cells and compute mean
cell_out = mean(cat(length(size(cell_in))+1, cell_in{:}), length(size(cell_in))+1);