% BLP estimation with demand side moments
load('data/100_3.mat')

% if dataset has zero shares, we need to replace the values
shares(shares == 0) = 10e-6;
prices(shares == 0) = 0;
prods(shares == 0, :) = 0;

% run the estimation procedure
[theta, fval] = blpdemand(prices, prods, shares)

