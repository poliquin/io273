% Create and save datasets
% ========================

% Problem 1.3
% -----------
rng(8675309);  % seed for reproducibility

% 10 markets with 3 products
[shares, prices, prods, profits, surplus, xi] = mktsim(3, 10);
save('data/10_3.mat', 'shares', 'prices', 'prods', 'profits', 'surplus', 'xi');
mkthist('tex/figs/hist_10_3.pdf', prices, profits, surplus, '(j,m)=(3,10)');

% 100 markets with 3 products
[shares, prices, prods, profits, surplus, xi] = mktsim(3, 100);
save('data/100_3.mat', 'shares', 'prices', 'prods', 'profits', 'surplus', 'xi');
mkthist('tex/figs/hist_100_3.pdf', prices, profits, surplus, '(j,m)=(3,100)');

% 100 markets with 5 products
[shares, prices, prods, profits, surplus, xi] = mktsim(5, 100);
save('data/100_5.mat', 'shares', 'prices', 'prods', 'profits', 'surplus', 'xi');
mkthist('tex/figs/hist_100_5.pdf', prices, profits, surplus, '(j,m)=(5,100)');

