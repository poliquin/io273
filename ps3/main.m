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
F = orderstat(X, 20, 160, 2);

% create a histogram of simulation values
f = figure('PaperPosition', [.1, .2, 6.2, 3.5], 'PaperSize', [6.4, 4]);
p1 = plot(20:2:160, mean(F, 2))
xlabel('Value')
ylabel('')
title('CDF of Value Distribution')
saveas(f, 'figs/ascending_values.pdf');

