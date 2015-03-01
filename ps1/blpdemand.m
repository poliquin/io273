function [theta, gammas, vcov, fval] = blpdemand(prices, prods, shares, ...
        cost, zc, prodcount, mktcount, usegrad, usecost, conduct)
    % BLP Estimation of model from BLP (1995)
    % Input arguments:
    %   prices = mXj by 1 vector of prices for each product-market combination
    %   prods  = mXj by 3 vector of product characteristics
    %   shares = mXj by 1 vector of market shares
    %   cost = j by 1 vector of firm-level marginal cost shifters
    %   zc   = j*m by 1 vector of market and product specific cost variables
    %   prodcount = number of products
    %   mktcount  = number of markets
    %   usegrad   = true means use analytic gradient of objective function in
    %               optimization routine
    %   usecost   = use cost-shifter moment condition in demand estimation
    %   conduct   = model of firm behavior \in {oligopoly, monopoly, perfect}
    % Outputs:
    %   theta = [alpha; beta; sigma_alpha]
    %   vcov = variance covariance matrix for theta
    %   fval = value of objective function evaluated at theta
    global deltas;

    %% construct matrix of BLP instruments for demand side
    % ------------------------------------------------------------------------
    Z = abs(eye(prodcount) - 1);   % matrix with 1 on off diagonal
    Z = repmat({Z}, mktcount, 1);  % selects other products in market
    Z = blkdiag(Z{:}) * prods;  % sum of other product characteristics
    % remove the first instrument, which does not vary across markets
    % because all markets have the same number of products.
    Z = Z(:, 2:3);

    % add marginal cost shifters to the demand-side instrument matrix
    cost = repmat(cost, mktcount, 1);
    if usecost
        Z = [cost, Z];
    end

    % this is the Nevo instrument
    totprice = repmat(eye(prodcount), 1, mktcount) * prices;
    avgprice = (1/(mktcount - 1))*(repmat(totprice, mktcount, 1) - prices);
    % uncomment below to add Nevo instrument
    %Z = [avgprice, Z];
    
    %% construct matrix of BLP instruments for the supply side
    % ------------------------------------------------------------------------
    K = [cost, zc];  % use the given cost-side variables
    cons = ones(prodcount * mktcount, 1);
    % the following would create BLP cost-side instruments, but our cost
    % variables are not endogenous...
    %K = abs(eye(prodcount) - 1);   % matrix with 1 on off diagonal
    %K = repmat({K}, mktcount, 1);  % selects other products in market
    %K = blkdiag(K{:});
    % supply-side instruments
    %K = [K * cost, K * zc];  % sum of other product cost shifters
    
    %% Set initial deltas to the solution of logit equation
    % ------------------------------------------------------------------------
    % find an initial value for delta using logit model
    share_mat = reshape(shares, prodcount, []);  % markets by products matrix
    out_share = 1 - sum(share_mat);  % share of outside option
    
    deltas = log(share_mat ./ repmat(out_share, 3, 1));
    deltas = reshape(deltas, 1, [])';
    % uncomment below to estimate the logit model
    %[coef, ~, ~] = logit(logit_shr, [prices, prods], [Z, prods]);
    
    %% draw 500 consumers for each market, held constant over simulation
    % ------------------------------------------------------------------------
    nu = lognrnd(0, 1, mktcount, 500);
    nu = kron(nu, ones(prodcount, 1));  % replicate draws for each product

    % initial weighting matrix
    % TODO: Change weighting matrix at each iteration
    W1 = inv([Z, prods]' * [Z, prods]);  % demand-side weights
    W2 = inv([cons, K]' * [cons, K]);    % supply-side weights
    W = blkdiag(W1, W2);                 % moment condition weights
    
    %% Run estimation routine
    % ------------------------------------------------------------------------
    tolerance = 1e-12;  % tolerance for inner loop, stricter than outer loop
    estimator = @(s) gmmobj(s, prods, Z, tolerance, conduct);
    options = optimset('Display', 'iter', 'TolFun', 1e-10);
    if usegrad  % use the gradient info in optimization routine
        options = optimset(options, 'GradObj', 'on');
        % uncomment below to check derivative against finite difference
        %options = optimset(options, 'DerivativeCheck', 'on');
        [s, fval, grad] = fminunc(estimator,  unifrnd(-1, 2), options);
    else
        [s, fval] = fminunc(estimator, lognrnd(0,1), options);
    end
    vcov = stderr(exp(s));
    
    % ------------------------------------------------------------------------
    function [fval, grad] = gmmobj(sigma, X, Z, innertol, conduct)
        % GMMOBJ Objective function for BLP random coefficients model
        % Input arguments:
        %   sigma = random coefficient on price, sigma_alpha
        %   X = matrix of product characteristics
        %   Z = matrix of BLP instruments
        %   innertol = tolerance for inner loop
        %   conduct = model of firm behavior {oligopoly, monopoly, perfect}
        % Outputs:
        %   fval = value of the objective function
        %   grad = gradient for speedy computation
       
        % sigma must be a positive number
        sigma = exp(sigma);

        %% Esimtate deltas
        % --------------------------------------------------------------------
        % create a simulator for market shares that can calculate shares for
        % different values of delta (the d variable); this is used to equate 
        % the observed shares with simulated shares and thereby find deltas.
        price_utility = sigma * bsxfun(@times, nu, prices);  % price disutility
        
        % reshape variables for inner loop
        price_utility_inner = reshape(price_utility, prodcount, []);
        deltas_inner = reshape(deltas, prodcount, []);
        shares_inner = reshape(shares, prodcount, []);
        sharefunc = @(d) deltashares(d, price_utility_inner, prodcount);
        
        % find deltas using the share simulator
        deltas_in = innerloop(deltas_inner, shares_inner, sharefunc, innertol);
        deltas = deltas_in(:);  % reshape deltas back
        
        % make sure deltas are defined, otherwise set high objective value
        if any(isnan(deltas)) == 1
            fval = 10e6;  % avoid these regions
            return
        end

        %% Estimate demand-side non-random coefficients and unobservables
        % --------------------------------------------------------------------
        [betas, xi] = ivreg(deltas, [prices, X], [Z, X], W1);
        
        %% Calculate markups and supply-side parameters
        % --------------------------------------------------------------------
        % there are 3 possible conduct models: perfect competition, oligopoly,
        % and perfect collusion (monopoly).
        if strcmp(conduct, 'monopoly')
            marks = collusion(sigma, nu, prices, deltas, shares, ...
                              prodcount, mktcount);
        elseif strcmp(conduct, 'perfect')
            marks = log(prices); 
        else  % oligopoly is the default
            marks = markup(sigma, nu, prices, deltas, shares, ...
                           prodcount, mktcount);
        end    
        [gammas, eta] = ivreg(marks, [cons, cost, zc], [cons, K], W2);
    
        %% Compute value of the objective function
        % --------------------------------------------------------------------
        % Calculate demand and supply-side moments
        dmom = [Z, X]' * xi;       % demand-side moment conditions
        smom = [cons, K]' * eta;   % supply-side moment conditions
        fval = [dmom; smom]' * W * [dmom; smom];  % stack 'em high
        
        if nargout > 1
            % find the jacobian, then calculate gradient of objective function
            jac = jacob(deltas, sigma, prices, nu, prodcount, mktcount);
            grad = 2*jac' * [Z, X, cons, K] * W * [dmom; smom];
        end
    
        % save latest parameter values
        theta = [betas; sigma];
    end
    % ------------------------------------------------------------------------
    function [vcov] = stderr(sigma)
        % STDERR Calculate standard errors for BLP parameter estimates
        %   This code follows Nevo's example code.

        % derivative of share delta function with respect to sigma
        D = jacob(deltas, sigma, prices, nu, prodcount, mktcount);
        % calculate first portion of matrix
        IVD = [Z, prods];  % full demand-side instrument matrix
        IVS = [cons, K];   % full supply-side instrument matrix
        QD = [prices, prods, D]' * IVD;
        QS = [cons, cost, zc]' * IVS;
        % calculate inverse of Gamma' * Gamma, see page 858 of BLP (1995)
        Q = blkdiag(QD, QS);
        B = inv(Q * W * Q');
        
        % find residuals for final parameters
        [~, xi] = ivreg(deltas, [prices, prods], [Z, prods], W1);
        marks = markup(sigma, nu, prices, deltas, shares, prodcount, mktcount);
        [~, eta] = ivreg(marks, [cons, cost, zc], [cons, K], W2);
        
        % calculate interaction of instrument matrix and residuals
        SD = IVD .* (xi * ones(1, size(IVD, 2)));
        SS = IVS .* (eta * ones(1, size(IVS, 2)));
        V1 = [SD, SS]' * [SD, SS];  % V1 from page 858 of BLP (1995)
        % covariance matrix using just V1
        vcov = B * Q * W * V1 * W * Q' * B;
    end
end
