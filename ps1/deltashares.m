function [shares] = deltashares(deltas, prutil, prodcount)
    % DELTASHARES Calculate market shares implied by given delta values
    %   Function f_j from BLP (1995), see equations 6.6 and 6.10
    % Input arguments:
    %   deltas = mXj vector of product specific components of utility function.
    %   prutil = Consumer specific disutility resulting from price, it is an
    %            j*m by N matrix with product/markets in rows and consumers
    %            in cols. 
    %   prodcount = Number of products per market.
    % Outputs:
    %   shares = Implied market shares, a mXj vector
   
%     % DYK: I think you are averaging before taking exponents;
%     % use each delta to calculate utility for all simulated consumers
%     u = bsxfun(@minus, deltas, prutil);
%     mean_u = mean(exp(u), 2);  % mean utility for each product
% 
%     % For each market, we take the mean over simulated consumers following
%     % equation 6.10, which is the simple simulator proposed by BLP (1995).
%     tops = reshape(mean_u, prodcount, []);  % products in rows, markets in cols
%     shares = bsxfun(@rdivide, tops, 1 + sum(tops));
% 
%     shares = reshape(shares, [], 1);  % stack the shares
    u = bsxfun(@minus, deltas, prutil);
    tops = exp(reshape(u,prodcount,[]));
    bottom = 1+sum(tops);
    pre_shares = reshape(bsxfun(@rdivide,tops,bottom),size(deltas,1),[]);
    shares = mean(pre_shares,2);
end

% prutil is the disutility each consumer gets from price of each product,
% it is a N by j*m matrix with consumers in rows and products in columns.
% calculate is as coefficients * prices' where coefficients are (1 + 1*nu).

