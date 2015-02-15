% Create and save datasets
% ========================

% Problem 1.3
% -----------
rng(8675309);  % seed for reproducibility

% 10 markets with 3 products
[shares, prices, products, profits] = mktsim(3, 10);
save('data/10_3.mat', 'shares', 'prices', 'products', 'profits');

% 100 markets with 3 products
[shares, prices, products, profits] = mktsim(3, 100);
save('data/100_3.mat', 'shares', 'prices', 'products', 'profits');

% 100 markets with 5 products
[shares, prices, products, profits] = mktsim(5, 100);
save('data/100_5.mat', 'shares', 'prices', 'products', 'profits');

