function [theta, fval] = blpdemand(prices, prods, shares)
    % BLP Estimation of model from BLP (1995)
    %
    
    prodcount = 3;  % TODO: replace with actual product count
    mktcount = 100; % TODO: replace with actual market count

    % find an initial value for delta using logit model
    share_mat = reshape(shares, prodcount, []);  % markets by products matrix
    out_share = 1 - sum(share_mat);  % share of outside option
    deltas = reshape(bsxfun(@minus, share_mat, out_share), [], 1);

    % construct Nevo instrument
    totprice = repmat(eye(prodcount, prodcount), 1, mktcount) * prices;
    avgprice = (1/(mktcount - 1))*(repmat(totprice, mktcount, 1) - prices);

    % construct matrix of BLP instruments
    Z = abs(eye(prodcount) - 1);   % matrix with 1 on off diagonal
    Z = repmat({Z}, mktcount, 1);
    Z = blkdiag(Z{:}) * prods;  % sum of other product characteristics
    Z = [Z];

    % draw consumers, held constant for all simulations
    nu = lognrnd(0, 1, 1500, 1);

    % initial weighting matrix
    W = ([prods, Z]' * [prods, Z]) \ eye(size([prods, Z], 2));

    options = optimset('Display', 'iter');
    estimator = @(s) gmmobj(s, deltas, prices, prods, Z, W, shares, nu);
    [sigma, fval] = fminunc(estimator, 1, options);

    function [fval, grad] = gmmobj(sigma, deltas, prices, X, Z, W, shares, nu)
        % GMMOBJ Objective function for BLP random coefficients model
        % Input arguments:
        %   theta = model parameters [beta; alpha; sigma_alpha]
        %   deltas = initial mXj by 1 vector of values for delta
        %   prices = mXj by 1 vector of prices
        %   X = matrix of product characteristics
        %   Z = matrix of BLP instruments
        %   shares = mXj by 1 vector of market shares
        % Outputs:
        %   fval = value of the objective function
        %   grad = gradient for speedy computation
       
        % create a simulator for market shares that can calculate shares for
        % different values of delta (the d variable); this is used to equate 
        % the observed shares with simulated shares and thereby find deltas.
        price_utility = sigma * nu * prices';  % consumer disutility from price
        sharefunc = @(d) deltashares(d, price_utility);  % share simulator
        
        % find deltas using the share simulator
        tolerance = 2e-8;
        deltas = innerloop(deltas, shares, sharefunc, tolerance);
        
        % estimate non-random coefficients and unobservables
        [betas, xi] = ivreg(deltas, [X, prices], [X, Z], W);
    
        % compute value of the objective function
        fval = xi' * [X, Z] * W * [X, Z]' * xi;
        
        if nargout > 2
            grad = 0;  % TODO: compute gradient of objective function
        end
    
        % save latest parameter values
        theta = [betas; sigma];
    end
end

