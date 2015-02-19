% Section 2, Question 1

load('data/100_3.mat');
% number of markets and products in each market
mkt_count = size(surplus, 1);
prod_count = size(prods, 1) / mkt_count;

% Unobserved and observed product characteristics
moment1 = mean(mean(bsxfun(@times, xi, prods), 2));
sprintf('E[xi_jm X_jm] = %0.5f', moment1)

% Unobserved characteristics and prices
moment2 = mean(xi .* prices);
sprintf('E[xi_jm p_jm] = %0.5f', moment2)

% Unobserved characteristics and mean prices in other markets
totprice = repmat(eye(prod_count, prod_count), 1, mkt_count) * prices;
avgprice = (1/(mkt_count - 1))*(repmat(totprice, mkt_count, 1) - prices);
moment3 = mean(xi .* avgprice);
sprintf('E[xi_jm bar{p}_jm] = %0.5f', moment3)

