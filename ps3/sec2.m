clear all

rng(8675309);
dir = mfilename('fullpath');
cd(dir(1:end-4));

%% Simulate dataset
N=100;  K=2; SIGMA=[1,1;1,2]; ALPHA=1; BETA=1;
[Y, X, P,A] = sim_dataset(N, K, SIGMA, ALPHA, BETA);

%% Run gibbs
[A,B,S] = gibbs(Y,[X,P],5000,[-10,-10],100);