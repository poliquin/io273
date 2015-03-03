function [shares] = deltashares(deltas, prutil, prodcount)
    % DELTASHARES Calculate market shares implied by given delta values
    %   Function f_j from BLP (1995), see equations 6.6 and 6.10
    % Input arguments:
    %   deltas = j by m matrix of product specific components of utility function.
    %   prutil = Consumer specific disutility resulting from price, it is a
    %            j by m*N matrix with products in rows and
    %            market/individuals in columns
    % Outputs:
    %   shares = Implied market shares, a mXj vector

    % simulate consumers following equation 6.10, which is the simple 
    % simulator proposed by BLP (1995).
    tops = exp(repmat(deltas, 1, 500) - prutil);
    % tops(tops==0) = 10e-20; % If deltas are zero, change it to small value
    bottom = 1 + sum(tops);
    shares = reshape(bsxfun(@rdivide, tops, bottom), prodcount ... 
        * size(deltas, 2), []);
end

% prutil is the disutility each consumer gets from price of each product,
% it is a N by j*m matrix with consumers in rows and products in columns.
% calculate it as coefficients * prices' where coefficients are sigma*nu.

