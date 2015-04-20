% Problem Set 3
% =============
% Do Yoon Kim, Chris Poliquin, David Zhang
% April 28, 2015

rng(8675309);
dir = mfilename('fullpath');
cd(dir(1:end-4));

%% Section 1, Question 1
% ----------------------------------------------------------------------------

X = dlmread('ascending_data.dat');
F = orderstat(X, 20, 160, 2, 0);

% create a histogram of simulation values
f = figure('PaperPosition', [.1, .2, 6.2, 3.5], 'PaperSize', [6.4, 4]);
p1 = plot(20:2:160, mean(F, 2))
xlabel('Value')
ylabel('')
title('CDF of Value Distribution')
saveas(f, 'figs/ascending_values.pdf');

%% Section 1, Question 2 (Haile & Tamer, 2003)
% ----------------------------------------------------------------------------

% upper bound is minimum of previously estimated values from order statistic
ub = min(F, [], 2);

% lower bound is maximum of estimated values, assuming min bid increment = 1
lb = max(orderstat(X, 20, 160, 2, 1), [], 2);

x = [20:2:160]';
f = figure('PaperPosition', [.1, .2, 6.2, 3.5], 'PaperSize', [6.4, 4]);
ax1 = axes('Parent', f, 'FontSize', 11);
area(x, lb, 'FaceColor', [0.5 0.9 0.6])
hold on
area(x, ub, 'FaceColor', [1.0 1.0 1.0])
%plot(x, mean(F, 2), 'r--')  % uncomment to also plot earlier cdf estimate
hold off
xlabel('Value', 'FontSize', 12);
ylabel('')
title('Bounds on CDF of Value Distribution', 'FontSize', 14)
saveas(f, 'figs/ascending_bounds.png');

%% Section 1, Question 3 (asymmetric private values, first price auction)
% ----------------------------------------------------------------------------

fpa = dlmread('fpa.dat');
[T M] = size(fpa);

FU = zeros(T, M);  % matrix to hold estimated values
bw = 12;  % bandwidth parameter for kernel function
% get psuedo private values for each bidder, i
for i = 1:M
    sprintf('Estimating values for bidder %i', i)
    bidi = fpa(:,i);          % player i bids
    noti = fpa;
    noti(:,i) = [];           % bidders other than player i
    maxj = max(noti, [], 2);  % max opponent bid
    
    % find value for player i in each auction
    FU(:,i) = bidshade(bidi, maxj, bw);
end

% use psuedo private values to estimate empirical cdf for each bidder
ui = mat2cell(FU, T, ones(1, M));
cdf = arrayfun(@(z) cell2mat(cellfun(@(x) mean(x < z), ui, 'UniformOutput', ...
                             false)), [0:2:160]', 'UniformOutput', false);
cdf = cell2mat(cdf);

% plot value distribution for each bidder
f = figure('PaperPosition', [.1, .2, 6.2, 3.5], 'PaperSize', [6.4, 4]);
p1 = plot(0:2:160, cdf);
xlabel('Value')
ylabel('')
title('CDF of Value Distributions, 4 Bidders')
saveas(f, 'figs/fpa_values.pdf');

