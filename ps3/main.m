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
lb = max(orderstat_n_2(X, 20, 160, 2, 1), [], 2);

x = [20:2:160]';
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
ui = mat2cell(FU, T, ones(1, M));
cdf = arrayfun(@(z) cell2mat(cellfun(@(x) mean(x < z), ui, 'UniformOutput', ...
                             false)), (0:2:250)', 'UniformOutput', false);
cdf = cell2mat(cdf);

% plot value distribution for each bidder
f = figure('PaperPosition', [.1, .2, 6.2, 3.5], 'PaperSize', [6.4, 4]);
p1 = plot(0:2:250, cdf);
xlabel('Value')
ylabel('')
title('CDF of Value Distributions, 4 Bidders')
saveas(f, 'figs/fpa_values.pdf');

% evaluate cdf of value distributions at 25th percentile of valuations
p25 = arrayfun(@(z) cell2mat(cellfun(@(x) mean(x < z), ui, 'UniformOutput', ...
                        false)), prctile(FU, 25)', 'UniformOutput', false);
p25 = cell2mat(p25)
% evaluate cdf of value distributions at 75th percentile of valuations
p75 = arrayfun(@(z) cell2mat(cellfun(@(x) mean(x < z), ui, 'UniformOutput', ...
                        false)), prctile(FU, 75)', 'UniformOutput', false);
p75 = cell2mat(p75)

%% Section 2, Question 1
% ----------------------------------------------------------------------------

% Simulate data
[Y, X, P, A] = sim_dataset(100, 2, [1, 1; 1, 2], 1, 1);

% Set location of priors
beta = 2;
alpha = 2;
sigma = [1,2,; 2,6];
% Set diffusion of priors
alpha_var = 1;
beta_var = 1;
v = 1;

for sim = 1:1
%% Step 1: Draw latent varaible w
% Draw latent variable w, which has two columns, from the truncated 
% multivariate normal distribution via rejection sampling
for i=1:100
    flag = 0; % initialize rejection flag
    while flag==0
        w(i,:) = mvnrnd(X(i,:)*beta(sim) + P(i,:)*alpha(sim), sigma(:,:,sim)); % draw from multivariate normal
        yhat = A(:,:,i)*w(i,:)'; % Find latent y
        % Implement indicator function for y
        yhat(yhat<0)=0; % Set negative latent y to 0
        yhat = all(yhat); % Set observed latent y to 1 if all elements are non-zero, or 0 otherwise
        % Reject if do not agree with data
        if Y(i)==yhat
            flag = 1;
        end
    end
end

% Transform to stacked form
X_stack = [X(:,1);X(:,2)];
P_stack = [P(:,1);P(:,2)];
w_stack = [w(:,1);w(:,2)];

X_reg = [X_stack,P_stack];
beta_reg = [beta(sim); alpha(sim)];
A_reg = [beta_var, 0; 0,alpha_var];

%% Step 2: Find posterior beta
beta_hat = inv(X_reg'*X_reg + A_reg)*(X_reg'*w_stack + A_reg*beta_reg);
beta(sim+1) = beta_hat(1,1);
alpha(sim+1) = beta_hat(2,1); 

%% Step 3: Find posterior sigma_hat
epsilon_hat =  w - (X*beta(sim+1) + P*alpha(sim+1));
for i =1:100
    ep_sum(:, :,i) = epsilon_hat(i,:)'*epsilon_hat(i,:);
end
ep_sum = sum(ep_sum,3);

inv_sigmahat = wishrnd(sigma(:,:,sim) + ep_sum, v + 100);
sigma(:, :, sim+1) = inv(inv_sigmahat);
end

