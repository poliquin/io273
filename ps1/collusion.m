function [marks] = collusion(sigma, nu, prices, deltas, shares, ...
        prodcount, mktcount)
    % COLLUSION  Calculate markups assuming perfect collusion.
    %   This function is used to calculate the dependent variable to estimate
    %   supply-side parameters assuming monopoly-like conduct.
    % Input arguments:
    %   sigma  = standard deviation of alpha
    %   nu     = consumers, j*m by N matrix where N is number of consumers
    %   prices = j*m by 1 vector of prices
    %   deltas = j*m by 1 vector of simulated deltas
    %   prodcount = number of products
    %   mktcount  = number of markets
    % Outputs:
    %   marks = j*m by 1 vector of markups
    
    N = size(nu, 2);  % number of simulated consumers
    
    %% calculate probability each consumer chooses each product
    % ------------------------------------------------------------------------
    % TODO: do not repeat code from deltashares.m, should share a function
    prutil = sigma * bsxfun(@times, nu, prices);  % price disutility
    u = bsxfun(@minus, deltas, prutil);
    tops = exp(reshape(u, prodcount, []));
    bottom = 1 + sum(tops);
    probs = reshape(bsxfun(@rdivide, tops, bottom), size(deltas, 1), []);
    
    %% (dsdp) partial derivative of share j with respect to price j
    % ------------------------------------------------------------------------
    dsdp = mean(-1 * nu .* probs .* (1-probs), 2);  % diagonal of delta matrix
 
    %% (dsdk) partial derivative of share j with respect to price k
    % ------------------------------------------------------------------------
    % group the observations into cells by market
    mktcells = mat2cell(probs, prodcount*ones(mktcount, 1), N);
    % get the possible combinations of 2 products (prodcount choose 2)
    combos = nchoosek(1:prodcount, 2)';
    % function that picks pairs of products from purchase probability matrix,
    % multiplies the shares and nu, and takes mean across consumers
    dsdk_mean = @(x, mc, mktnu) mean(prod([mc(combos(:,x), :); mktnu]));
    dsdk = @(mc, mktnu) arrayfun(@(x) dsdk_mean(x, mc, mktnu), 1:prodcount);
   
    %% calculate markups
    % ------------------------------------------------------------------------
    % delta matrix element (i, k) is partial derivative of share i with respect
    % to price k; it is a diagonal matrix in our oligopoly case.
    marks = zeros(mktcount*prodcount, 1);  % place to store markups
    dsize = [prodcount, prodcount];  % dimensions of delta matrix within market
    for mkt=1:mktcount
        dmatrix = zeros(prodcount);  % delta matrix for this market
        % first row corresponding to this market
        i = prodcount*(mkt - 1) + 1;
        % all rows corresponding to this market
        mrows = i:(i + prodcount - 1);
        % find dsdk and put elements on upper triangle
        consumers = nu(i, :);
        dsdk_mkt = dsdk(mktcells{mkt}, consumers);
        dmatrix(sub2ind(dsize, combos(1,:), combos(2,:))) = dsdk_mkt;
        % convert triangular matrix to symmetric matrix
        dmatrix = dmatrix + triu(dmatrix)';
        % put elements from dsdp on diagonal
        dmatrix(logical(eye(prodcount))) = dsdp(mrows);
        % calculate the markup
        marks(mrows) = prices(mrows) - dmatrix\shares(mrows);
    end
    marks = log(marks);  % equation 3.6 from BLP (1995), page 854
end
