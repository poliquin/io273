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
beta    = .99;
RC      = 10;
theta1  = .05;
theta30 = .3;
theta31 = .5;
theta32 = .2;

EV = ev(beta, RC, theta1, theta30, theta31, theta32);

%% Section 2, Question 5
% ----------------------------------------------------------------------------
[xt, it] = sim_data(EV,beta, RC, theta1, theta30, theta31);


%% Section 3, Question 1
% ----------------------------------------------------------------------------
%Step 1: estimate transition probabilities
[theta31hat,theta31hatse,theta32hat,theta32hatse] = transitionprob(xt,it);

%Step 2: Use nested fixed point to find theta1
options = optimset('Display', 'iter', 'TolFun', 10e-10);
[theta1hat,~,~,~,~,hess] = fminunc(@(theta1hat) -1 * rust(theta1hat, theta31hat, theta32hat, RC, beta, xt, it), rand, options);
se = sqrt(diag(inv(hess)));
sprintf('3.1\ntheta1 = %f (%f), \ntheta31 = %f (%f), \ntheta32 = %f(%f)', theta1hat, se,theta31hat,theta31hatse, theta32hat, theta32hatse)


