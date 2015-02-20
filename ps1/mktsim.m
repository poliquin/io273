function [shares, prices, products, profits, surplus, xi, W] = mktsim(j, m)
    %MKTSIM Draw j products in m markets for BLP simulation.
    %   Simulate 500 individuals per market for m markets with j products,
    %   creating characteristics, cost shifters, prices, and market shares.
    % Input arguments:
    %   j = number of products
    %   m = number of markets
    % Outputs:
    %   jXm by 1 vector of observed (simulated) market shares
    %   jXm by 1 vector of prices
    %   jXm by 3 matrix of product characteristics
    %   jXm by 1 vector of profits
    %     m by 1 vector of total consumer surplus for each market
    %   jXm by 1 vector of unobserved characteristics
    %     j by 1 vector of firm-specific marginal cost shifters
    % Population parameters
    N = 2000;          % number of consumers
    ALPHA = 1;         % population constant price sensitivity
    SIGMA = 1;         % variance of price sensitivity
    BETA = [5; 1; 1];  % characteristic coefficients
    GAMMA = [2; 1; 1]; % marginal cost coefficients
    options = optimoptions('fsolve', 'Display', 'iter');
    options = optimoptions(options, 'MaxIter', 800);
    
    % Pre-allocate space for resulting dataset
    shares = zeros(j * m, 1);    % market shares, simulated below
    prices = zeros(j * m, 1);    % prices, solved below
    products = zeros(j * m, 3);  % characteristics, drawn below
    profits = zeros(j * m, 1);   % profit for each product/market, found below
    surplus = zeros(m, 1);       % total consumer surplus, calculated below
    xi = zeros(j * m, 1);        % unobserved characteristics

    % Marginal cost shifter for firm j across all markets
    W = normrnd(0, 1, j, 1);

    for k=1:m  % simulate each market individually
        
        % Consumer tastes
        nu = lognrnd(0, 1, N, 1);     % consumer specific price sensitivity
        epsilon = evrnd(0, 1, N, j);  % type 1 extreme value
        alpha_i = ALPHA + SIGMA * nu; % random coefficient on price
       
        while 1
            % Model does not always converge, so we need to re-run the
            % simulation until it converges in each market. We can re-draw
            % the characteristics and costs to help model find workable values.
            
            % Product characteristics: constant, uniform, and standard normal
            X = [ones(j, 1), unifrnd(0, 1, j, 1), normrnd(0, 1, j, 1)];
            xi_k = normrnd(0, 1, j, 1);   % unobserved characteristics

            % Market-specific marginal cost variables
            Z = normrnd(0, 1, j, 1);
            eta = normrnd(0, 1, j, 1);
            MC = 0.2 * exp([ones(j, 1), W, Z] * GAMMA + eta); 

            % Nash equilibrium in market k
            firmobj = @(P) equilibrium(P, BETA, X, MC, alpha_i, xi_k);
            P0 = ones(j, 1);
            [P, fval] = fsolve(firmobj, P0, options);
            
            % Calculate utilities (consumers in rows, products in columns)
            ploss = bsxfun(@times, alpha_i, P');
            U = repmat(BETA'*X', N, 1) - ploss + repmat(xi_k', N, 1) + epsilon;

            [simulated_shares, simulated_surplus] = simshare(U);
            % Check if model solved with positive prices and shares
            if fval < 10^-4 & P > 0 & sum(simulated_shares) > 0
                % found a workable solution for market k
                % solve for profits
                profit_k = (P - MC) .* (N * simulated_shares);
                break
            end
        end

        % Save simulation results for market k to dataset
        shares(1+(k-1)*j : k*j) = simulated_shares;
        prices(1+(k-1)*j : k*j) = P;
        products(1+(k-1)*j : k*j, :) = X;
        profits(1+(k-1)*j : k*j) = profit_k;
        surplus(k) = simulated_surplus;
        xi(1+(k-1)*j : k*j) = xi_k;
    end  % end simulation for market k
end

function inner_share = link_fun(P, BETA, X, alpha_i, xi)
    % LINK_FUN First part of market share equation (integrating over epsilon)
    %   Solves the first part of eq 6.7 in BLP, 1995.
    % Outputs:
    %   Function for computing first part of market share equation.
    disutil = bsxfun(@times, alpha_i, P');
    posutil = BETA' * X' + xi';
    denom = mean(exp(bsxfun(@minus, posutil, disutil)));
    inner_share = denom / (1 + sum(denom));
end

function foc = equilibrium(P, BETA, X, MC, alpha_i, xi)
    % EQUILIBRIUM Solving firm's first order condition.
    %   Computes first-order condition for Nash equilibrium in single market.
    % Outputs:
    %   Value of the first-order condition.
    share = link_fun(P, BETA, X, alpha_i, xi);
    % Derivatives of market share with respect to prices (eq 6.9a in BLP, 1995)
    inner_ds = bsxfun(@times, -1*alpha_i, share.*(1 - share));
    ds = mean(inner_ds);
    
    foc = share' + (P - MC) .* ds';
end

function [shares, surplus] = simshare(U)
    %MKTSHARE Calculate market share of each product.
    %   Finds optimal choice for each consumer, then aggregates over consumers
    %   to produce simulation market shares. This differs from firm's problem
    %   because firms don't see epsilon when setting prices.
    % Input Arguments:
    %   U = N by j matrix of consumer utilities
    % Outputs:
    %   Vector of market shares for each product.
    
    % Find product that gives maximum utility
    [maxU, maxI] = max(U, [], 2);
    maxI(logical(maxU < 0)) = 0;  % these consumers choose outside option
    maxU(logical(maxU < 0)) = 0;  % utility from outside option is zero

    % Tabulate market shares from consumer choices
    [n, product_count] = size(U);
    shares = zeros(product_count, 1);
    for j=1:product_count
        shares(j, 1) = sum(maxI == j) / n;
    end

    % Sum consumer surplus across all consumers
    surplus = sum(maxU);
end

