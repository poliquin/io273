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
z(shares == 0) = 0;

%% Run the estimation procedure for demand and supply side
% ----------------------------------------------------------------------------
runs = 10;
% run the procedure multiple times so that we get multiple starting values
estimates = zeros(runs, 9);  % 9 columns: fval, 5 demand & 3 supply coefs
% TODO: add variance-covariance matrix for gamma
variances = zeros(runs, 64); % 64 columns to hold re-shaped 8x8 matrix
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
vcov = reshape(variances(minidx, :), 8, 8);
vcov_theta = vcov(1:5, 1:5);
vcov_gamma = vcov(6:8, 6:8);

% calculate bias of the estimates
bias_theta = theta - [-1; 5; 1; 1; 1];
% we multiplied MC by 0.2 in simulation, so the simulation constant term in
% the marginal cost function is actually log(0.2) + 2
bias_gamma = gammas - [log(0.2) + 2; 1; 1];

% print the best estimates
disp(' ')
disp('(3, 100) Demand Side Estimates')
disp(table(num2str(theta, '%3.3f &'), ...
           num2str(sqrt(diag(vcov_theta)), '(%3.3f) &'), ...
           num2str(bias_theta, '%4.3f \\\\'), ...
            'VariableNames', {'Coef' 'SE' 'Bias'}, ...
            'RowNames', {'Price  &', 'X1  &', 'X2  &', 'X3  &', 'Sigma  &'}))
fprintf('Objective function: %i', fval)
disp(' ')

disp('(3, 100) Supply Side Estimates')
disp(table(num2str(gammas, '%3.3f &'), ...
           num2str(sqrt(diag(vcov_gamma)), '(%3.3f) &'), ...
           num2str(bias_gamma, '%4.3f \\\\'), ...
            'VariableNames', {'Coef' 'SE' 'Bias'}, ...
            'RowNames', {'Cons &', 'W  &', 'Z  &'}))
fprintf('Objective function: %i', fval)
disp(' ')
