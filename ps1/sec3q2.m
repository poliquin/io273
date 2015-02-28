% BLP estimation of supply and demand side

%% Load the (3, 100) dataset
% ----------------------------------------------------------------------------
load('data/100_3.mat')
prodcount = 3;
mktcount = 100;

% if dataset has zero shares, we need to replace the values
tabulate(shares == 0)
shares(shares == 0) = 10e-6;
prices(shares == 0) = 0;
prods(shares == 0, :) = 0;
cost(shares == 0, :) = 0;

%% Run the estimation procedure for demand and supply side
% ----------------------------------------------------------------------------
runs = 1;
% run the procedure multiple times so that we get multiple starting values
estimates = zeros(runs, 9);  % 9 columns: fval, 5 demand & 3 supply coefs
% TODO: add variance-covariance matrix for gamma
variances = zeros(runs, 25); % 25 columns to hold re-shaped 5x5 matrix
for i=1:runs  % run multiple times
    fprintf('Run %1.0f of %1.0f', i, runs)
    [theta, gammas, vcov, fval] = blpdemand(prices, prods, shares, ...
            cost, z, prodcount, mktcount, true, true);
    estimates(i, :) = [fval, theta', gammas'];
    variances(i, :) = reshape(vcov, 1, []);
    disp(' ')
end

% show all of the estimates
disp('Full results')
disp(estimates)
disp('Variance of the estimates across runs')
disp(var(estimates, 1))

% find best set of coefficients and standard errors (smallest fval)
[fval, minidx] = min(estimates(:, 1));
theta  = estimates(minidx, 2:6)';
gammas = estimates(minidx, 7:9)';
vcov = reshape(variances(minidx, :), 5, 5);

% calculate bias of the estimates
bias_theta = theta - [-1; 5; 1; 1; 1];
% TODO: find real bias, we multiplied MC by 0.2 in simulation!
bias_gamma = gammas - [2; 1; 1];

% print the best estimates
disp(' ')
disp('(3, 100) Demand Side Estimates')
disp(table(num2str(theta, '%3.3f &'), ...
           num2str(sqrt(diag(vcov)), '(%3.3f) &'), ...
           num2str(bias_theta, '%4.3f \\\\'), ...
            'VariableNames', {'Coef' 'SE' 'Bias'}, ...
            'RowNames', {'Price  &', 'X1  &', 'X2  &', 'X3  &', 'Sigma  &'}))
fprintf('Objective function: %i', fval)
disp(' ')

% TODO: include standard errors for supply-side estimates
disp('(3, 100) Supply Side Estimates')
disp(table(num2str(gammas, '%3.3f &'), ...
           num2str(bias_gamma, '%4.3f \\\\'), ...
            'VariableNames', {'Coef' 'Bias'}, ...
            'RowNames', {'Cons &', 'W  &', 'Z  &'}))
fprintf('Objective function: %i', fval)
disp(' ')


