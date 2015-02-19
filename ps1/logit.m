function [coef, predict, resid] = logit(Y, X, Z, W)
    % LOGIT Basic logit model
    %   Estimate logit discrete model with instrumental variables regression.
    % Input arguments:
    %   Y = dependent variable vector, shares minus outside share
    %   X = design matrix of [price | product characteristics]
    %   Z = matrix of instruments, [instruments | product characteristics]
    %   W = optional weight matrix

    % calculate coefficients using instrumental variables regression
    if nargin < 4
        [coef resid] = ivreg(Y, X, Z);
    else
        [coef resid] = ivreg(Y, X, Z, W);
    end
    % find predicted values
    predict = X * coef;
end    

