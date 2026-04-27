function vectorSum = bbar_vecsum(data)
% COMPUTEVECTORSUM Computes the vector sum (magnitude) for each row of XYZ accelerometer data
%
% INPUT:
%   data  - Nx3 matrix where columns are [X, Y, Z]
%
% OUTPUT:
%   vectorSum - Nx1 vector of magnitudes for each sample
%
% USAGE:
%   vectorSum = bbar_vecsum(data);

    % Validate input
    if size(data, 2) ~= 3
        error('Input must be an Nx3 matrix with columns [X, Y, Z]');
    end

    vectorSum = sqrt(sum(data.^2, 2));

end