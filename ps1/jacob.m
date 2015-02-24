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
    prutil = sigma * bsxfun(@times, nu, prices);  % price disutility
    utils = exp(bsxfun(@minus, deltas, prutil));  % exponentiated utility
    
    % group the observations into cells by market and simulated consumer
    % so that purchase probability can be calculated within consumer
    utils = mat2cell(utils, prodcount*ones(mktcount, 1), ones(N, 1));
    share = @(x) x ./ (1 + sum(x));  % calculate purchase probability
    probs = cellfun(share, utils, 'UniformOutput', false);
    probs = cell2mat(probs);
    
    %% (dsdj) partial derivative of share j with respect to delta j
    % ------------------------------------------------------------------------
    dsdj = mean(probs .* (1-probs), 2);  % diagonal of derivative w.r.t. delta

    %% (dsdm) partial derivative of share j with respect to delta m
    % ------------------------------------------------------------------------
    % group the observations into cells by market and simulated consumer
    consumers = mat2cell(probs, prodcount*ones(mktcount, 1), ones(N, 1));
    % multiply each purchase probability by the other purchase probabilities
    cprods = cellfun(@(x) kron(x, x), consumers, 'UniformOutput', false);
    % calculate the mean of these products across consumers
    mus = mean(cell2mat(cprods), 2);
    % group the means by market, products in row and a single column
    mus = mat2cell(mus, prodcount*prodcount*ones(mktcount, 1), 1);
    % reshape the means into a product by product matrix in which
    % element (j, m) equals partial derivative of share j with respect
    % to delta m for j not equal to m.
    mu_reshaper = @(x) reshape(x, prodcount, prodcount);
    dsdm = cellfun(mu_reshaper, mus, 'UniformOutput', false);

    %% (dsdk) partial derivative of share j with respect to sigma k
    % ------------------------------------------------------------------------
    % derivative entails multiplying product prices by shares, within consumer
    prices_x_shares = bsxfun(@times, prices, probs);
    prices_x_shares = mat2cell(prices_x_shares, ...
                               prodcount*ones(mktcount, 1), ones(N, 1));
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
    % group dsdj and dsdk elements into cells by market
    dsdj = mat2cell(dsdj, prodcount*ones(mktcount, 1), 1);
    dsdk = mat2cell(dsdk, prodcount*ones(mktcount, 1), 1);
    for mkt=1:mktcount
        % for each market, put elements of cell in dsdj on diagonal of dsdm
        dsdm{mkt, 1}(logical(eye(prodcount))) = dsdj{mkt, 1};
        % gradient calculation for this market
        dsdm{mkt, 1} = -1 * inv(dsdm{mkt, 1}) * dsdk{mkt, 1};
    end
    % multiply negative inverse of above by dsdk
    jac = cell2mat(dsdm);

end
