% Problem Set 2
% =============
% Do Yoon Kim, Chris Poliquin, David Zhang
% March 30, 2015

rng(8675309);
addpath('derivest')  % files from John D'Errico's DERIVEST suite

%% 2.2(1) Generate and save data for the entry game

[mrkts, costs, firms, entry] = sim_markets(3, 100);
save('data/entry.mat', 'mrkts', 'costs', 'firms', 'entry');
% create a histogram of simulation values
f = figure('PaperPosition', [.1, .2, 6.2, 3.5], 'PaperSize', [6.4, 4]);
subplot(1,3,1)
p1 = histogram(mrkts(:,2));
xlim([0 3])
title('Number of entrants')
hold on
subplot(1,3,2)
p2 = histogram(mrkts(:,3:end));
title('Realized profits')
subplot(1,3,3)
p3 = scatter(mrkts(:,2), mean(mrkts(:,3:end), 2));
xlabel('Number of entrants')
ylabel('Average profits')
title('Profits and Entrants')
saveas(f, 'figs/sim.pdf');


%% 2.2(2) Maximum likelihood estimation with correct model

% draw standard normals for the simulation estimator
[M, F] = size(firms);  % number of markets and potential entrants
draws = normrnd(0, 1, 100, M*F);

theta = [1, 1, 1];  % true, known alpha, beta, delta
options = optimset('Display', 'iter', 'TolFun', 10e-10);

% likelihood function with mu = x(1) and sigma = x(2)
like = @(x, ord) berry(mrkts, firms, entry, x(1), x(2), theta, draws, ord);
initial = [unifrnd(-1, 4), unifrnd(0, 3)];
sprintf('Starting estimation at mu = %f and sigma = %f', initial)

% estimation assuming entry in order of profitability 2.2(2a)
[x1, fval1] = fminsearch(@(x) -1 * like(x, 'ascend'), initial, options);
[hess1, ~] = hessian(@(x) -1 * like(x, 'ascend'), x1);
se1 = sqrt(diag(inv(hess1)));
sprintf('2.2(2a)\nmu = %f (%f)\nsigma = %f (%f)', x1(1), se1(1), x1(2), se1(2))

% estimation assuming entry in reverse order 2.2(2b)
[x2, fval2] = fminsearch(@(x) -1 * like(x, 'descend'), initial, options);
[hess2, ~] = hessian(@(x) -1 * like(x, 'descend'), x2);
se2 = sqrt(diag(inv(hess2)));
sprintf('2.2(2b)\nmu = %f (%f)\nsigma = %f (%f)', x2(1), se2(1), x2(2), se2(2))


%% 2.3 Estimate mean costs of entry using moment inequality estimator
NumSims = 100; % Number of simulations
theta = [1, 1, 1, 1];  % true, known alpha, beta, delta, sigma
% draw u
for i=1:NumSims
    u(:,:,i) = normrnd(0, theta(4), size(mrkts, 1), size(mrkts, 2) - 2);
end
% Find muhat
[muhat] = fminsearch(@(mu) moment_inequalities(theta, mu, mrkts, firms, entry, u), ...
                     unifrnd(-1, 4), options);

