function signal = bbar_zigzagScan(block)
    % Perform zigzag scanning of a 2D block to create 1D signal
    [M, N] = size(block);
    signal = zeros(M*N, 1);
    
    index = 1;
    for sum = 0:(M+N-2)
        if mod(sum,2) == 0  % going up
            for i = min(sum,M-1):-1:max(0,sum-N+1)
                j = sum - i;
                signal(index) = block(i+1,j+1);
                index = index + 1;
            end
        else  % going down
            for i = max(0,sum-N+1):min(sum,M-1)
                j = sum - i;
                signal(index) = block(i+1,j+1);
                index = index + 1;
            end
        end
    end
end