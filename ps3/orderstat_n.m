function [F] = orderstat(X, small, large, stepsize, increment)
    % ORDERSTAT Estimate CDF of private values, F(v) from ascending auction
    %           data using distribution for n-1 order statistic, G(v).
    % Input arguments:
    %   X = T by 2 matrix with auction rows and #bidders, price in columns
    %   small = lowest value of v to estimate F(v) for
    %   large = highest value of v to estimate F(v) for
    %   stepsize = size of increments between small and large
    %   increment = minimum bid increment
    % Outputs:
    %   F = values of F(v) for auctions with n bids

    syms x;  % distribution of n-1 order statistic with n bidders
    ostat = @(n, g, x) x^n - g;  % g = G(v)
    
    unq = unique(X(:,1));  % number of bidders that participate in auctions
    nu = size(unq, 1);
    grp = arrayfun(@(i) X(X(:,1)==i, 2), unq, 'UniformOutput', false);
    
    steps = ceil((large - small) / stepsize);
    F = zeros(steps, nu);
    
    i = 0;
    for v = small:stepsize:large
        i = i + 1;
        % get cdf of n-1 order statistic at value v for auctions with N bidders
        g = cellfun(@(bids) mean(bids + increment < v), grp);

        % find roots for n-1 order statistic, which equal F(v)
        % this function solves expression for n-1 order statistic for F(v),
        % given a number of bidders n and a value for G(v).
        f = arrayfun(@(t) fzero(@(x) ostat(unq(t), g(t), x), [0 1]), 1:nu);
        F(i,:) = f;
    end

    F = F(1:i,:);
end
