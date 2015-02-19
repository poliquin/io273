function [shares] = deltashares(deltas, prutil, prodcount)
    % DELTASHARES Calculate market shares implied by given delta values
    %   Function f_j from BLP (1995), see equations 6.6 and 6.10
    % Input arguments:
    %   deltas = mXj vector of product specific components of utility function.
    %   prutil = Consumer specific disutility resulting from price, it is an
    %            N by j*m matrix with consumers in rows and products in cols. 
    %   prodcount = Number of products per market.
    % Outputs:
    %   shares = Implied market shares, a mXj vector
   
    % use each delta to calculate utility for all simulated consumers
    u = bsxfun(@minus, deltas, prutil);
    mean_u = mean(exp(u), 2);  % mean utility for each product

    % For each market, we take the mean over simulated consumers following
    % equation 6.10, which is the simple simulator proposed by BLP (1995).
    tops = reshape(mean_u, prodcount, []);  % products in rows, markets in cols
    shares = bsxfun(@rdivide, tops, 1 + sum(tops));

    shares = reshape(shares, [], 1);  % stack the shares
end

% prutil is the disutility each consumer gets from price of each product,
% it is a N by j*m matrix with consumers in rows and products in columns.
% calculate is as coefficients * prices' where coefficients are (1 + 1*nu).

