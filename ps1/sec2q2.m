% BLP estimation with demand side moments

runs = 1;  % number of times to run BLP algorithm

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

%% Run the estimation procedure using demand side moments
% ----------------------------------------------------------------------------
% run the procedure multiple times so that we get multiple starting values
estimates = zeros(runs, 6);  % 6 columns: fval plus coefficients
variances = zeros(runs, 25); % 25 columns to hold re-shaped 5x5 matrix
for i=1:runs  % run six times
    fprintf('Run %1.0f of %1.0f', i, runs)
    [theta, ~, vcov, fval] = blpdemand(prices, prods, shares, cost, z, ...
            prodcount, mktcount, true, false, false, 'oligopoly');
    estimates(i, :) = [fval, theta'];
    variances(i, :) = reshape(vcov, 1, []);
    disp(' ')
end

% show all of the estimates
disp('Full results')
disp(estimates)

% find best set of coefficients and standard errors (smallest fval)
[fval, minidx] = min(estimates(:, 1));
theta = estimates(minidx, 2:end);
vcov  = reshape(variances(minidx, :), 5, 5);

% calculate bias of the estimates
bias = theta' - [-1; 5; 1; 1; 1];

% print the best estimates
disp(' ')
disp('(3, 100) Coefficients and Standard Errors (Demand-Side Moments)')
disp(table(num2str(theta', '%3.3f &'), num2str(sqrt(diag(vcov)), '(%3.3f) &'), ...
            num2str(bias, '%4.3f   \\\\'), ...
            'VariableNames', {'Coef' 'SE' 'Bias'}, ...
            'RowNames', {'Price  &', 'X1  &', 'X2  &', 'X3  &', 'Sigma  &'}))
fprintf('Objective function: %i', fval)
disp(' ')


%% Load the (3, 10) dataset
% ----------------------------------------------------------------------------
clear all
load('data/10_3.mat')
prodcount = 3;
mktcount = 10;

% if dataset has zero shares, we need to replace the values
tabulate(shares == 0)
shares(shares == 0) = 10e-6;
prices(shares == 0) = 0;
prods(shares == 0, :) = 0;

%% Run the estimation procedure using demand side moments
% ----------------------------------------------------------------------------
% run the procedure multiple times so that we get multiple starting values
estimates = zeros(runs, 6);  % 6 columns: fval plus coefficients
variances = zeros(runs, 25); % 25 columns to hold re-shaped 5x5 matrix
for i=1:runs  % run six times
    fprintf('Run %1.0f of %1.0f', i, runs)
    [theta, ~, vcov, fval] = blpdemand(prices, prods, shares, cost, z, ...
            prodcount, mktcount, true, false, false, 'oligopoly');
    estimates(i, :) = [fval, theta'];
    variances(i, :) = reshape(vcov, 1, []);
    disp(' ')
end

% show all of the estimates
disp('Full results')
disp(estimates)

% find best set of coefficients and standard errors (smallest fval)
[fval, minidx] = min(estimates(:, 1));
theta = estimates(minidx, 2:end);
vcov  = reshape(variances(minidx, :), 5, 5);

% calculate bias of the estimates
bias = theta' - [-1; 5; 1; 1; 1];

% print the best estimates
disp(' ')
disp('(3, 10) Coefficients and Standard Errors (Demand-Side Moments)')
disp(table(num2str(theta', '%3.3f &'), num2str(sqrt(diag(vcov)), '(%3.3f) &'), ...
            num2str(bias, '%4.3f \\\\'), ...
            'VariableNames', {'Coef' 'SE' 'Bias'}, ...
            'RowNames', {'Price  &', 'X1  &', 'X2  &', 'X3  &', 'Sigma  &'}))
fprintf('Objective function: %i', fval)
disp(' ')
