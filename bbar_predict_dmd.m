% Dynamic Mode Decomposition
%
% [X_hat_t, time_evolution] = bbar_predict_dmd(phi, delta, b0, t)
%
% After computing the Phi and Delta using bbar_dmd(), future state
% predictions can be obtained by:
%
% X_hat_t = Phi * Delta^t * b0        , b0 are the initial conditions.
%
% 
% Input:
% Phi   = Dynamic modes as complex numbers
% Delta = Diagonal of eigenvalues
% b0    = Initial conditions
%
% by Zakaria Djebbara, Feb. 2025
function [X_hat_t, time_evolution] = bbar_predict_dmd(phi, D, b0, t)
    
time_evolution = zeros(length(b0), length(t));   % Initializing timeseries

for j = 1:length(b0)
    time_evolution(j, :) = b0(j) * exp(log(D(j,j)) * t);  % Mode's time dynamics
end

X_hat_t = phi * time_evolution;

end