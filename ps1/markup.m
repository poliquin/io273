function [marks] = markup(sigma, nu, prices, deltas, shares, prodcount, mktcount)
    % MARKUP Calculate markups from oligopoly pricing equation.
    %   This function is used to calculate the dependent variable to estimate
    %   supply-side parameters. See step 3 on page 863 of BLP (1995).
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
    
    %% calculate markups
    % ------------------------------------------------------------------------
    % delta matrix element (i, k) is partial derivative of share i with respect
    % to price k; it is a diagonal matrix in our oligopoly case.
    marks = zeros(mktcount*prodcount, 1);  % place to store markups
    dsize = [prodcount prodcount];  % dimensions of delta matrix within market
    for mkt=1:mktcount
        dmatrix = zeros(prodcount);  % delta matrix for this market
        % first row corresponding to this market
        i = prodcount*(mkt - 1) + 1;
        % all rows corresponding to this market
        mrows = i:(i + prodcount - 1);
        % put elements from dsdp on diagonal
        dmatrix(logical(eye(prodcount))) = dsdp(mrows);
        % calculate the markup
        marks(mrows) = prices(mrows) - dmatrix\shares(mrows);
    end
    marks = log(marks);  % equation 3.6 from BLP (1995), page 854
end
