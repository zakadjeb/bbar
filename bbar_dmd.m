% Dynamic Mode Decomposition
%
% [phi, omega, lambda, b, Xdmd, timeDyn, a, E] = bbar_dmd(x1, x2, n)
%
% DMD decomposes the data into spatial and temporal modes. These are
% eigenvectors and eigenvalues, respectively. It's mainly based on the SVD
% decomposition and follows this structure:
%
% Xt+1 = A Xt           , Assuming linear structure in data
% 
% Xt+1 = X' = U S V*    , Basically SVD.
%
% U* X' = U* A U S V*   , Left-multiply U-transposed, and as V* V = I
%
% U* X' V = U* A U S    , S is just a diagonal matrix
%
% U* X' V S^-1 = U* A U , This gives an approximation of A:
%
% U* A U = Ã            , We know that A left-multiplied with the eigenvectors
% gives us the eigenvectors left-multiplied with the diagonal eigenvalues:
%
% Ã W = W D             , W: the eigenmatrix. D(delta): the eigenvalues.
%
% Phi = X' V S^-1 W     , that's the DMD.
% 
% INPUTS:
% x1        = x(1:end-1). M x N, M: temporal evolution. N: columnvector of data
% x2        = x(2:end). M x N, M: temporal evolution. N: columnvector of data
% n         = number of modes to compute.
% dt        = Time between frames
%
% OUTPUT:
% phi       = Dynamic Modes as complex numbers
% omega     = Continuous time-series for modes
% lambda    = Discrete time-series for modes
% b         = Vector of magnitude of modes
% Xdmd      = Reconstructed X , phi * b0^time
% timeDyn   = mode activations over time
% a         = data-driven mode weights over time (pinv(phi)*x1)
% E         = model measure (approximation error using Frobenius norm)
%
% by Zakaria Djebbara, Feb. 2025
% 
% Updates:
% 04-01-26: Added a model measure for selecting number of modes.
%

function [phi, omega, lambda, b, Xdmd, timeDyn, a, E] = bbar_dmd(x1, x2, n, dt)

if nargin < 4; dt = 1; end

% Computing the SVD
[U, S, V] = svd(x1,"econ"); % Computing U, S, V

% Computing the variables
U = U(:,1:n);                 % Reducing the number of U-modes
V = V(:,1:n);                 % Reducing the number of V-value
S = S(1:n,1:n);               % Reducing the number of scaling-values
S_inverse = diag(1./diag(S)); % Computing the inverse Sigma

% Computing the A_tilde (an approximation of A)
A_tilde = conj(U') * x2 * V * S_inverse; % Computing the Ã. 
[W,D] = eig(A_tilde);         % Computing the eigenmatrix and eigenvalues.

% Computing the DMD modes
phi = x2 * V * S_inverse * W; % Computing the DMD
lambda = diag(D);             % The discrete eigenvalues (activation)
omega = log(lambda)/dt;       % The continuous eigenvalues

% Computing the magnitude of each mode
b = pinv(phi) * x1(:,1);

% DMD reconstruction
m = size(x1,2);
t = (0:m-1)*dt;
timeDyn = b .* exp(omega * t);
Xdmd = phi * timeDyn;

% Data-driven mode weights over time (projection coefficients)
a = pinv(phi) * x1;

% Computing a model-measure for selecting no. of modes.
num = norm(x1 - Xdmd, 'fro')^2;
den = norm(x1, 'fro')^2;
E  = 1 - (num / den);
end