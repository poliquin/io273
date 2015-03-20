% Problem Set 2
% =============
% Do Yoon Kim, Chris Poliquin, David Zhang
% March 2015

rng(8675309);

%% Generate and save data for the entry game
[markets, costs, firms] = sim(3, 100);
save('data/entry.mat', 'markets', 'costs', 'firms');
% create a histogram of simulation values
f = figure('PaperPosition', [.1, .2, 6.2, 3.5], 'PaperSize', [6.4, 4]);
subplot(1,3,1)
p1 = histogram(markets(:,2));
xlim([0 3])
title('Number of entrants')
hold on
subplot(1,3,2)
p2 = histogram(markets(:,3:end));
title('Realized profits')
subplot(1,3,3)
p3 = scatter(markets(:,2), mean(markets(:,3:end), 2));
xlabel('Number of entrants')
ylabel('Average profits')
title('Profits and Entrants')
saveas(f, 'figs/sim.pdf');


