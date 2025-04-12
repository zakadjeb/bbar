function [RP, RQA] = bbar_recurrenceAnalysis(sig, varargin)
    % RECURRENCEANALYSIS Compute recurrence plot and RQA metrics for a signal
    %
    % Inputs:
    %   signal: Input time series (vector)
    %   params: Structure with fields:
    %       - embeddingDim: Embedding dimension (default: 2)
    %       - delay: Time delay/space offset for embedding (default: 1)
    %       - threshold: Distance threshold for recurrence (default: 0.1)
    %       - norm: Distance norm ('euclidean' or 'maximum', default: 'euclidean')
    %
    % Outputs:
    %   RP: Recurrence plot matrix
    %   RQA: Structure containing RQA metrics
    
    % Set default parameters if not provided
    p = inputParser;
    addRequired(p, 'sig', @isnumeric);
    addOptional(p, 'embeddingDim', 2, @isnumeric);
    addOptional(p, 'delay', 1, @isnumeric);
    addOptional(p, 'maxThreshold', 0.1, @isnumeric);
    % addOptional(p, 'minThreshold', eps+eps, @isnumeric);
    addOptional(p, 'norm', 'euclidean', @ischar);
    addOptional(p, 'minLineLength', 2, @isnumeric);
    parse(p,sig,varargin{:});
    
    % Phase space reconstruction using time/space-delay embedding
    N = length(sig);
    M = N - (p.Results.embeddingDim-1)*p.Results.delay;
    Y = zeros(M, p.Results.embeddingDim);
    if M > 3000; error('Matrix is too big. Consider using bbar_avgWindow().'); end

    for i = 1:p.Results.embeddingDim
        Y(:,i) = sig((1:M) + (i-1)*p.Results.delay);
    end
    
    % Compute distance matrix
    RP = sparse(zeros(M));
    for i = 1:M
        for j = 1:M
            if strcmp(p.Results.norm, 'euclidean')
                dist = sqrt(sum((Y(i,:) - Y(j,:)).^2));
            elseif strcmp(p.Results.norm, 'maximum')  % maximum norm
                dist = max(abs(Y(i,:) - Y(j,:)));
            end
            % RP(i,j) = (p.Results.minThreshold <= dist) & dist <= p.Results.maxThreshold;
            RP(i,j) = dist <= p.Results.maxThreshold;
        end
    end
    
    % Compute RQA metrics
    RQA = computeRQAMetrics(RP,p.Results.minLineLength);
end

function RQA = computeRQAMetrics(RP, minLineLength)
    % Compute various RQA metrics from the recurrence plot
    
    N = size(RP, 1);
    
    % Recurrence Rate (RR)
    RQA.RR = sum(sum(RP)) / (N^2);
    
    % Determinism (DET): Percentage of recurrence points forming diagonal lines
    [diagonalLines, totalDiagPoints] = findDiagonalLines(RP, minLineLength);
    if totalDiagPoints > 0
        RQA.DET = totalDiagPoints / sum(sum(RP));
    else
        RQA.DET = 0;
    end
    
    % Average diagonal line length (L)
    if ~isempty(diagonalLines)
        RQA.L = mean(diagonalLines);
    else
        RQA.L = 0;
    end
    
    % Maximum diagonal line length (Lmax)
    if ~isempty(diagonalLines)
        if ~isempty(max(diagonalLines(diagonalLines < size(RP,1))))
            RQA.Lmax = max(diagonalLines(diagonalLines < size(RP,1)));
            % RQA.Lmax = max(diagonalLines);
        else
            RQA.Lmax = 0;
        end
    else
        RQA.Lmax = 0;
    end
    
    % Entropy of diagonal line lengths (ENTR)
    if ~isempty(diagonalLines)
        [uniqueLengths, ~, counts] = unique(diagonalLines);
        p = counts / sum(counts);
        RQA.ENTR = -sum(p .* log(p));
    else
        RQA.ENTR = 0;
    end
    
    % Laminarity (LAM): Percentage of recurrence points forming vertical lines
    [verticalLines, totalVertPoints] = findVerticalLines(RP, minLineLength);
    if totalVertPoints > 0
        RQA.LAM = totalVertPoints / sum(sum(RP));
    else
        RQA.LAM = 0;
    end
    
    % Trapping time (TT): Average vertical line length
    if ~isempty(verticalLines)
        RQA.TT = mean(verticalLines);
    else
        RQA.TT = 0;
    end
end

function [lines, totalPoints] = findDiagonalLines(RP, minLength)
    % Find diagonal lines in the recurrence plot
    N = size(RP, 1);
    lines = [];
    totalPoints = 0;
    
    % Check each diagonal
    for k = -(N-minLength):(N-minLength)
        diagonal = diag(RP, k);
        lineStarts = find(diff([0; diagonal; 0]) == 1);
        lineEnds = find(diff([0; diagonal; 0]) == -1) - 1;
        lengths = lineEnds - lineStarts + 1;
        
        % Keep only lines longer than minLength
        validLines = lengths >= minLength;
        lines = [lines; lengths(validLines)];
        totalPoints = totalPoints + sum(lengths(validLines));
    end
end

function [lines, totalPoints] = findVerticalLines(RP, minLength)
    % Find vertical lines in the recurrence plot
    N = size(RP, 1);
    lines = [];
    totalPoints = 0;
    
    % Check each column
    for j = 1:N
        column = RP(:,j);
        lineStarts = find(diff([0; column; 0]) == 1);
        lineEnds = find(diff([0; column; 0]) == -1) - 1;
        lengths = lineEnds - lineStarts + 1;
        
        % Keep only lines longer than minLength
        validLines = lengths >= minLength;
        lines = [lines; lengths(validLines)];
        totalPoints = totalPoints + sum(lengths(validLines));
    end
end