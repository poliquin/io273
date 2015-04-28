% Problem Set 3
% =============
% Do Yoon Kim, Chris Poliquin, David Zhang
% April 28, 2015
clear all

rng(8675309);
dir = mfilename('fullpath');
cd(dir(1:end-4));

%% Section 1, Question 1
% ----------------------------------------------------------------------------

X = dlmread('ascending_data.dat');
F = orderstat(X, 20, 160, 2, 0);

f = figure('PaperPosition', [.1, .2, 6.2, 3.5], 'PaperSize', [6.4, 4]);
p1 = plot(20:2:160, mean(F, 2))
xlabel('Value')
ylabel('')
title('CDF of Value Distribution')
saveas(f, 'figs/ascending_values.pdf');

%% Section 1, Question 2 (Haile & Tamer, 2003)
% ----------------------------------------------------------------------------
% Set smoothing parameter, which increases with N for consistency
rho = -1/600*size(X,1);
% Assume first price English auction
% upper bound is min of estimated values of the nth order statistic
y = orderstat_n(X,20,160,1,0);
numerator = y.*exp(y*rho);
denominator = sum(exp(y*rho),2);
ub = sum(numerator./repmat(denominator,1,3), 2);

% lower bound is maximum of estimated values of the n-1th order statistic
y = orderstat(X,20,160,1,1);
numerator = y.*exp(y*-rho);
denominator = sum(exp(y*-rho),2);
lb = sum(numerator./repmat(denominator,1,3), 2);

x = [20:1:160]';
f = figure('PaperPosition', [.1, .2, 6.2, 3.5], 'PaperSize', [6.4, 4]);
ax1 = axes('Parent', f, 'FontSize', 11);
area(x, ub, 'FaceColor', [0.5 0.9 0.6])
hold on
area(x, lb, 'FaceColor', [1.0 1.0 1.0])
%plot(x, mean(F, 2), 'r--')  % uncomment to also plot earlier cdf estimate
hold off
xlabel('Value', 'FontSize', 12);
ylabel('')
title('Bounds on CDF of Value Distribution', 'FontSize', 14)
saveas(f, 'figs/ascending_bounds.png');

%% Section 1, Question 3 (asymmetric private values, first price auction)
% ----------------------------------------------------------------------------

fpa = dlmread('fpa.dat');
[T, M] = size(fpa);

FU = zeros(T, M);  % matrix to hold estimated values
bw = 12.8;  % bandwidth parameter for kernel function
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
uir = mat2cell(FU, ones(1, T), 4);
cdf = @(z) mean(cell2mat(arrayfun(@(x) all(cell2mat(x) < z), ...
                                  uir, 'UniformOutput', false)));
% evaluate cdf at 16 combinations of 25th or 75th percentiles of values
tiles = {[25 75], [25 75], [25 75], [25 75]};
[w, x, y, z] = ndgrid(tiles{:});
combos = [w(:) x(:) y(:) z(:)];
numcmb = size(combos, 1);

cdfvalues = zeros(numcmb, 9);
for i = 1:numcmb
   quin = diag(prctile(FU, combos(i,:)))';  % get desired quintile of marginal
   cdfvalues(i,:) = [combos(i,:) quin, cdf(quin)];
end
sprintf('CDF Under Asymmetric Private Values')
cdfvalues

% Symmetric, independent, private values
v = cell2mat(arrayfun(@(x) symmetry(x, fpa, 12.8), fpa(:), ...
                      'UniformOutput', false));
pv = arrayfun(@(x) mean(v <= x), 0:2:150, 'UniformOutput', false);
f = figure('PaperPosition', [.1, .2, 6.2, 3.5], 'PaperSize', [6.4, 4]);
p1 = plot(0:2:150, cell2mat(pv))
xlabel('Value')
ylabel('')
title('CDF of Symmetric Value Distribution')
saveas(f, 'figs/symmetric_values.pdf');

%% Section 2, Question 1
% ----------------------------------------------------------------------------
