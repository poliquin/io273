function [shares] = deltashares(deltas, prutil, prodcount)
    % DELTASHARES Calculate market shares implied by given delta values
    %   Function f_j from BLP (1995), see equations 6.6 and 6.10
    % Input arguments:
    %   deltas = m by j matrix of product specific components of utility function.
    %   prutil = Consumer specific disutility resulting from price, it is a
    %            j by m*N matrix with products in rows and
    %            market/individuals in columns
    %   prodcount = Number of products per market.
    % Outputs:
    %   shares = m*j by N matrix of shares by individuals
    
    N_customers = size(prutil, 2) / size(deltas, 2); % Number of consumers
    J_markets = size(deltas, 2);
    
    % This function returns the conditional on nu market shares (BLP 6.6)
    % To obtain the shares for each market, you need to take the mean over 
    % simulated consumers following equation 6.10, which is the simple 
    % simulator proposed by BLP (1995).
    tops = exp(repmat(deltas,1,N_customers) - prutil);
    tops(tops==0) = 10e-5;
    bottom = 1 + sum(tops);
    shares = reshape(bsxfun(@rdivide, tops, bottom), prodcount ... 
        * J_markets, []);
    %shares = reshape(mean(pre_shares,2),prodcount,[]);
end

% prutil is the disutility each consumer gets from price of each product,
% it is a N by j*m matrix with consumers in rows and products in columns.
% calculate it as coefficients * prices' where coefficients are sigma*nu.

