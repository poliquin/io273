% Problem Set 4
% =============
% Do Yoon Kim, Chris Poliquin, David Zhang
% May 15, 2015
clear all

rng(8675309);
dir = mfilename('fullpath');
cd(dir(1:end-4));

%% Section 2, Question 3
% ----------------------------------------------------------------------------
% Define true parameters
beta    = .99;
RC      = 10;
theta1  = .05;
theta30 = .3;
theta31 = .5;
theta32 = .2;
% Calculate expected continuation value
EV = ev(beta, RC, theta1, theta30, theta31, theta32);

%% Section 2, Question 4
% Plot continuation value
EVplot = EV(1:11,:);
EVplotx = linspace(0,10,11);
f = figure('PaperPosition', [.1, .2, 6.2, 3.5], 'PaperSize', [6.4, 4]);
plot(EVplotx,EVplot(:,1)','--',EVplotx,EVplot(:,2),':')
title('Plot of Choice-specific Value Functions')
xlabel('x')
ylabel('EV')
legend('i = 0', 'i = 1')
saveas(f, 'figs/EV.pdf');


%% Section 2, Question 5
% ----------------------------------------------------------------------------
% Simulate data
[xt, it] = sim_data(EV,beta, RC, theta1, theta30, theta31);
% Generate summary statistics and plots
xtit = zeros(10000,2);
count = 1;
for i=1:size(xt,1)
    for j=1:size(xt,2)
        xtit(count, 1) = xt(i,j);
        xtit(count, 2) = it(i,j);
        count = count + 1;
    end
end
histogram(xtit(xtit(:,2)==1,1))
title('Histogram of x when i=0')
xlabel('x')
ylabel('freq')
saveas(f, 'figs/histx.pdf');
sprintf('2.5\nMin = %f, \n25th Pctile = %f, \nMean = %f, \nMedian = %f, \n75th Pctile = %f, \nMax=%f',...
    min(xtit(xtit(:,2)==1,1)), prctile(xtit(xtit(:,2)==1,1),25),mean(xtit(xtit(:,2)==1,1)),median(xtit(xtit(:,2)==1,1)), prctile(xtit(xtit(:,2)==1,1),75), max(xtit(xtit(:,2)==1,1)))


%% Section 3, Question 1
% ----------------------------------------------------------------------------
%Step 1: estimate transition probabilities
[theta31hat,theta31hatse,theta32hat,theta32hatse] = transitionprob(xt,it);

%Step 2: Use nested fixed point to find theta1
options = optimset('Display', 'iter', 'TolFun', 10e-10);
[theta1hat,~,~,~,~,hess] = fminunc(@(theta1hat) -1 * rust(theta1hat, theta31hat, theta32hat, RC, beta, xt, it), rand, options);
se = sqrt(diag(inv(hess)));
sprintf('3.1\ntheta1 = %f (%f), \ntheta31 = %f (%f), \ntheta32 = %f(%f)', theta1hat, se,theta31hat,theta31hatse, theta32hat, theta32hatse)


%% Section 3, Question 3
% ----------------------------------------------------------------------------
% columnize choices and transitions
renew = reshape(it(1:999, :), [], 1);
trans = reshape(diff(xt), [], 1);
% estimate transition probabilities
theta3 = tabulate(trans(renew == 0));
theta3 = theta3(:, 3);





