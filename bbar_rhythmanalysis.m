function [subRP, subRQA] = bbar_rhythmanalysis(img, N, M, varargin)
    % imgSUBDIVISIONRP Subdivide img and compute recurrence plots
    %
    % Inputs:
    %   img: Input img (2D matrix)
    %   N: Number of subdivisions in rows
    %   M: Number of subdivisions in columns
    %   params: Structure with RP parameters (passed to recurrenceAnalysis)
    %
    % Outputs:
    %   subRP: Cell array of recurrence plots for each subdivision
    %   subRQA: Cell array of RQA metrics for each subdivision
    
    p = inputParser;
    addRequired(p, 'img', @ismatrix);
    addRequired(p, 'N', @isnumeric);
    addRequired(p, 'M', @isnumeric);
    addOptional(p, 'gamma', 1/5, @isnumeric);
    addOptional(p, 'signalizing', 'zigzag', @ischar);

    addOptional(p, 'embeddingDim', 2, @isnumeric);
    addOptional(p, 'delay', 1, @isnumeric);
    addOptional(p, 'maxThreshold', 0.1, @isnumeric);
    % addOptional(p, 'minThreshold', eps+eps, @isnumeric);
    addOptional(p, 'norm', 'euclidean', @ischar);
    addOptional(p, 'minLineLength', 2, @isnumeric);

    parse(p,img,N,M,varargin{:});

    % Check if img is RGB and convert to grayscale if needed
    if size(img, 3) == 3
        img = rgb2gray(img)/255;
    end

    % Add white noise to avoid high correlation between subsequent zeros
    noise = remapnumbers(randn(size(img)),0,1);
    img = img + noise;
    img = img/max(max(img));

    if p.Results.gamma > 1; gamma = gamma / 1; end
    
    % Resize image according to gamma
    img = imresize(img,p.Results.gamma);
    [height, width] = size(img);
    
    % Calculate subdivision dimensions
    subHeight = floor(height/N);
    subWidth = floor(width/M);
    
    % Initialize cell arrays for results
    subRP = cell(N, M);
    subRQA = cell(N, M);
    
    % Process each subdivision
    for i = 1:N
        for j = 1:M
            % Displaying progress...
            warning('Now computing subdivision N: %d, M:%d.',i,j);

            % Calculate subdivision boundaries
            rowStart = (i-1)*subHeight + 1;
            rowEnd = min(i*subHeight, height);
            colStart = (j-1)*subWidth + 1;
            colEnd = min(j*subWidth, width);
            
            % Extract subdivision
            subimg = img(rowStart:rowEnd, colStart:colEnd);
            
            % Convert subdivision to 1D signal
            if strcmp(p.Results.signalizing, 'reshape')
                % Method 1
                signal = reshape(subimg,1,[]);
            elseif strcmp(p.Results.signalizing, 'row')
                % Method 2: Row-wise mean
                signal = mean(subimg, 2);
            elseif strcmp(p.Results.signalizing, 'column')
                % Method 3: Column-wise mean
                signal = mean(subimg, 1)';
            elseif strcmp(p.Results.signalizing, 'zigzag')
                % Method 4: Zigzag scanning
                signal = bbar_zigzagScan(subimg);
            end
            
            % Compute recurrence plot and RQA metrics
            [rp, rqa] = bbar_recurrenceAnalysis(signal,...
                'embeddingDim',p.Results.embeddingDim,...
                'delay',p.Results.delay,...
                'maxThreshold',p.Results.maxThreshold,...
                'norm',p.Results.norm,...
                'minLineLength',p.Results.minLineLength);

            % Store results
            subRP{i,j} = rp;
            subRQA{i,j} = rqa;
        end
    end
    warning('Done.');
end



