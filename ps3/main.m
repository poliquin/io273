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
F = orderstat(X, 34, 100, 2);
plot(34:2:100, mean(F, 2))

