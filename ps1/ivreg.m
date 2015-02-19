function [betas, xi] = ivreg(deltas, X, Z, W)
    % IVREG instrumental variables regression for linear parameters.
    %   See page 5 of the appendix to Nevo's practitioners' guide.
    % Input arguments:
    %   deltas = mXj by 1 vector of deltas computed from contraction mapping
    %   X = design matrix of terms common to all individuals, i.e. from delta
    %   Z = instrument matrix
    %   W = weight matrix
    % Outputs:
    %   betas = linear model parameters, same for all people, i.e. from delta
    %   xi = unobservables implied by estimate of deltas
    
    % user did not supply weight matrix, so create one
    if nargin < 4
        W = (Z' * Z) \ eye(size(Z, 2));
    end

    betas = (X' * Z * W * Z' * X) \ (X' * Z * W * Z' * deltas);
    xi = deltas - X * betas;
end

