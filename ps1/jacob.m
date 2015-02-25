function [jac] = jacob(deltas, sigma, prices, nu, prodcount, mktcount)
    % JACOB  Compute jacobian of market share inversion
    %    The derivative of the share invesion (which produces the deltas) is
    %    required to calculate the gradient of the GMM objective function.
    %    This formula comes from Nevo's practitioner guide; it requires the
    %    calculation of three separate derivatives:
    %      (1) of share j with respect to delta j
    %      (2) of share j with respect to delta m
    %      (3) of share j with respect to sigma (of random coefficient, nu)
    % Input arguments:
    %    deltas = m*j by 1 vector of deltas from the contraction mapping
    %    sigma = standard deviation of random coefficient on price
    %    prices = m*j by 1 vector of prices for each product in each market
    %    nu = m*j by N matrix with consumers in columns, products in rows
    %    prodcount = number of products in each market
    %    mktcount = number of markets
    % Outputs:
    %    jac = derivative of share inversion w.r.t non-linear parameter

    N = size(nu, 2);  % number of simulated consumers
    
    %% calculate probability each consumer chooses each product
    % ------------------------------------------------------------------------
    % TODO: do not repeat code from deltashares.m, should share a function
    prutil = sigma * bsxfun(@times, nu, prices);  % price disutility
    u = bsxfun(@minus, deltas, prutil);
    tops = exp(reshape(u, prodcount, []));
    bottom = 1 + sum(tops);
    probs = reshape(bsxfun(@rdivide, tops, bottom), size(deltas, 1), []);
    
    %% (dsdj) partial derivative of share j with respect to delta j
    % ------------------------------------------------------------------------
    dsdj = mean(probs .* (1-probs), 2);  % diagonal of derivative w.r.t. delta

    %% (dsdm) partial derivative of share j with respect to delta m
    % ------------------------------------------------------------------------
    % group the observations into cells by market
    mktcells = mat2cell(probs, prodcount*ones(mktcount, 1), N);
    % get the possible combinations of 2 products (prodcount choose 2)
    combos = nchoosek(1:prodcount, 2)';
    % function that picks pairs of products from purchase probability matrix,
    % multiplies the shares, and takes mean across consumers
    dsdm_mean = @(x, mc) mean(prod(mc(combos(:,x), :)));
    dsdm_product = @(mc) arrayfun(@(x) dsdm_mean(x, mc), 1:prodcount);
    dsdm = cellfun(dsdm_product, mktcells, 'UniformOutput', false); 
    
    %% (dsdk) partial derivative of share j with respect to sigma k
    % ------------------------------------------------------------------------
    % derivative entails multiplying product prices by shares, within consumer
    prices_x_shares = bsxfun(@times, prices, probs);
    prices_x_shares = mat2cell(prices_x_shares, prodcount*ones(mktcount, 1), N);
    % sum the shares multiplied by prices across products within each consumer
    share_summer = @(x) repmat(sum(x), prodcount, 1);
    csums = cellfun(share_summer, prices_x_shares, 'UniformOutput', false);
    % subtract sums from prices as specified in derivative formula
    csums = bsxfun(@minus, prices, cell2mat(csums));
    % multiply nu, shares, and the above element-wise (i.e. for each consumer)
    dsdk = nu .* probs .* csums;
    % take the mean across consumers to get desired derivative, a vector
    dsdk = mean(dsdk, 2);

    %% full derivative matrix
    % ------------------------------------------------------------------------
    jac = zeros(mktcount*prodcount, 1);
    dsize = [prodcount prodcount];
    for mkt=1:mktcount
        dmatrix = zeros(prodcount);
        % first row corresponding to this market
        i = prodcount*(mkt - 1) + 1;
        % all rows corresponding to this market
        mrows = i:(i + prodcount - 1);
        % place elements from dsdm on upper triangle
        dmatrix(sub2ind(dsize, combos(1,:), combos(2,:))) = dsdm{mkt};
        % convert triangular matrix to symmetric matrix
        dmatrix = dmatrix + triu(dmatrix)';
        % put elements from dsdj on diagonal
        dmatrix(logical(eye(prodcount))) = dsdj(mrows);
        % gradient calculation for this market
        jac(mrows) = -1 * dmatrix\dsdk(mrows);
    end 

end
