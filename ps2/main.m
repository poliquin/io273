% Problem Set 2
% =============
% Do Yoon Kim, Chris Poliquin, David Zhang
% March 30, 2015

rng(8675309);

%% 2.2(1) Generate and save data for the entry game

[mrkts, costs, firms, entry] = sim(3, 100);
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
draw = normrnd(0, 1, 100, M*F);

% likelihood function with mu = x(1) and sigma = exp(x(2))
theta = [1, 1, 1];  % true, known alpha, beta, delta
like = @(x) berry(mrkts, firms, entry, x(1), exp(x(2)), theta, draw);

options = optimset('Display', 'iter', 'TolFun', 10e-10);
[x, fval] = fminsearch(@(x) -1 * like(x), unifrnd(-1, 4, 2, 1), options);

sprintf('mu = %f\nsigma = %f', x(1), exp(x(2)))

%% 2.3 Estimate mean costs of entry using moment inequality estimator
NumSims = 100;
theta = [1, 1, 1,1];  % true, known alpha, beta, delta, sigma
% draw u
for i=1:NumSims
    u(:,:,i) = normrnd(0, theta(4), size(mrkts, 1), size(mrkts,2)-2);
end
% Find muhat
options = optimset('Display', 'iter', 'TolFun', 10e-10);
[muhat] = fminsearch(@(mu) moment_inequalities(theta, mu, mrkts, firms, u), 4, options);

