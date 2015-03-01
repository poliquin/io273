% Compute marginal costs using most preferred estimates
% ----------------------------------------------------------------------------
% assumes the preferred estimates are saved in data/supply_side_est.mat

load('data/100_3.mat')
load('data/supply_side_est.mat')
% cost-side variables
K = [ones(300, 1), repmat(cost, 100, 1), z];

% calculate true marginal costs
mc_true = 0.2 * exp(K * [2; 1; 1] + eta);

% calculate estimated marginal costs
mc_hat = 0.2 * exp(K * gammas);

% plot estimation error for products with non-zero share
mc_error = mc_hat - mc_true;

f = figure('PaperPosition', [.1, .2, 6, 3.5], 'PaperSize', [6.2, 4]);
h1 = histogram(mc_error(shares ~= 0), 25);  % 25 bins
title('Difference between estimated and true marginal cost')
annotation(f, 'textbox', [0.1,0.01,0.5,0.04], ...
           'String', {'(J, M) = (3, 100)'}, ...
           'LineStyle', 'none', ...
           'FontSize', 10);
saveas(f, 'tex/figs/hist_mc.pdf');

