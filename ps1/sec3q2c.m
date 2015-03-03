function [] = sec3q2c()
% SEC3Q2C  Estimate demand/supply model under different conduct assumptions

runs = 1;  % number of times to run each model

data = load('data/100_3.mat');
cost = data.cost;
eta = data.eta;
prices = data.prices;
prods = data.prods;
%profits = data.profits;
shares = data.shares;
%surplus = data.surplus;
%xi = data.xi;
z = data.z;

prodcount = 3;
mktcount = 100;

% if dataset has zero shares, we need to replace the values
tabulate(shares == 0)
shares(shares == 0) = 10e-6;
prices(shares == 0) = 0;
prods(shares == 0, :) = 0;
cost(shares == 0, :) = 0;
% cost-side variables
K = [ones(300, 1), repmat(cost, 100, 1), z];

%% Perfect collusion model
% ----------------------------------------------------------------------------
[theta, gammas, se_theta, se_gamma, collusion_mc] = runblp(runs, 'monopoly');
collusion = [theta, se_theta; gammas, se_gamma];

%% Perfect competition model
% ----------------------------------------------------------------------------
[theta, gammas, se_theta, se_gamma, competition_mc] = runblp(runs, 'perfect');
competition = [theta, se_theta; gammas, se_gamma];

%% Oligopoly model
% ----------------------------------------------------------------------------
[theta, gammas, se_theta, se_gamma, oligopoly_mc] = runblp(runs, 'oligopoly');
oligopoly = [theta, se_theta; gammas, se_gamma];


%% Create a table with all the estimates
% ----------------------------------------------------------------------------
disp(' ')
disp('Estimates under different conduct assumptions')
disp(table(num2str(oligopoly(:, 1), '%3.3f &'), ...
           num2str(oligopoly(:, 2), '(%3.3f) &'), ...
           num2str(collusion(:, 1), '%3.3f &'), ...
           num2str(collusion(:, 2), '(%3.3f) &'), ...    
           num2str(competition(:, 1), '%3.3f &'), ...
           num2str(competition(:, 2), '(%3.3f) \\'), ...
         'VariableNames', {'OCoef' 'OSE' 'MCoef', 'MSE', 'CCoef', 'CSE'}, ...
         'RowNames', {'Price  &', 'X1 &', 'X2 &', 'X3 &', 'Sigma  &', ...
                      'Cons  &', 'W &', 'Z &'} ...
    ))
disp(' ')

%% Create a plot of the marginal cost distributions 
% ----------------------------------------------------------------------------
f = figure('PaperPosition', [.1, .2, 6, 3.5], 'PaperSize', [6.2, 4]);
ksdensity(oligopoly_mc(:,2));
hold on
ksdensity(collusion_mc(:,2));
ksdensity(competition_mc(:,2));
legend('Oligopoly', 'Collusion', 'Perfect Competition')
title('Marginal costs under different conduct assumptions')
annotation(f, 'textbox', [0.1,0.01,0.5,0.04], ...
           'String', {'(J, M) = (3, 100)'}, ...
           'LineStyle', 'none', ...
           'FontSize', 10);
saveas(f, 'tex/figs/conduct_comparison.pdf');


%% Run BLP Algorithm
% ----------------------------------------------------------------------------
function [theta, gammas, se_theta, se_gamma, mc] = runblp(runs, conduct)
    % RUNBLP  Run BLP algorithm assuming given conduct.
    %
    %
    disp('---------------------------------------')
    fprintf('Estimating %i runs of %s model\n', runs, conduct)
    disp('---------------------------------------')
    % run the procedure multiple times so that we get multiple starting values
    estimates = zeros(runs, 9);  % 9 columns: fval, 5 demand & 3 supply coefs
    variances = zeros(runs, 64); % 64 columns to hold re-shaped 8x8 matrix
    for i=1:runs  % run multiple times
        fprintf('Run %1.0f of %1.0f', i, runs)
        [theta, gammas, vcov, ~, ~, fval] = blpdemand(prices, prods, shares, ...
                cost, z, prodcount, mktcount, true, true, true, conduct);
        estimates(i, :) = [fval, theta', gammas'];
        variances(i, :) = reshape(vcov, 1, []);
        disp(' ')
    end
    
    % find best set of coefficients and standard errors (smallest fval)
    [~, minidx] = min(estimates(:, 1));
    theta  = estimates(minidx, 2:6)';
    gammas = estimates(minidx, 7:9)';
    vcov = reshape(variances(minidx, :), 8, 8);
    
    % standard errors
    se_theta = sqrt(diag(vcov(1:5, 1:5)));
    se_gamma = sqrt(diag(vcov(6:8, 6:8)));
    
    % calculate true and estimated marginal costs
    mc_true = 0.2 * exp(K * [2; 1; 1] + eta);
    mc_hat = 0.2 * exp(K * gammas);
    mc_error = mc_hat - mc_true;
    mc = [mc_true, mc_hat, mc_error];
end

end
