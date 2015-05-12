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