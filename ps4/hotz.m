function [theta] = hotz(xt, it, RC, theta3, beta, opts)
    % HOTZ  Estimate theta_1 using Hotz & Miller (1996) inversion
    % Input arguments
    %   xt = matrix of states (mileage of bus in rows, buses in columns)
    %   it = renewal decisions (0 is continue, 1 is renew, buses in columns)
    %   RC = scalar engine replacement cost
    %   theta3 = 3 by 1 vector with estimated transition probabilities
    %   beta = scalar discount rate
    % Output arguments:
    %   theta = scalar estimate of theta_1

    % estimate conditional choice probabilities
    prob = ccp(xt, it);
    % prob(1) is the probability of replacing when xt = 0, and
    % prob(2) is the probability of replacing when xt = 1, and so on...
    
    x = xt(:);  % stack the observations
    % calculate difference in log probabilities (continue - renew)
    pdiff = log(1 - prob(x + 1)) - log(prob(x + 1));
    
    % calculate difference in values (continue - renew) at theta1
    cons = -beta*(theta3(1)*log(prob(x + 1)) + ...
                  theta3(2)*log(prob(x + 2)) + ...
                  theta3(3)*log(prob(x + 3)) ...
                 ) + RC + beta*log(prob(1));
         
    vdiff = @(theta) cons - (theta * x);
    
    % estimate theta using minimum distance
    theta = fminsearch(@(t) norm(pdiff - vdiff(t)), unifrnd(0, 1), opts);
end
