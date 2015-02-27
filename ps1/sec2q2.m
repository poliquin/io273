% BLP estimation with demand side moments
load('data/100_3.mat')
prodcount = 3;
mktcount = 100;

% if dataset has zero shares, we need to replace the values
tabulate(shares == 0)
shares(shares == 0) = 10e-6;
prices(shares == 0) = 0;
prods(shares == 0, :) = 0;

%% run the estimation procedure
% ----------------------------------------------------------------------------
% last parameter is whether to use jacobian during minimization of the
% objective function; doing so reduces time to convergence.
[theta, vcov, fval] = blpdemand(prices, prods, shares, cost, ...
        prodcount, mktcount, true);
disp('Coefficients:')
theta
disp('Standard Errors:')
sqrt(diag(vcov))

