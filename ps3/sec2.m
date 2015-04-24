clear all

rng(8675309);
dir = mfilename('fullpath');
cd(dir(1:end-4));

%% Simulate dataset
N=100;  K=2; SIGMA=[1,1;1,2]; ALPHA=1; BETA=1;
[Y, X, P,As] = sim_dataset(N, K, SIGMA, ALPHA, BETA);

%% Run gibbs using three different starting points and 10000 runs
runs = 3000; inits = 1;
alpha = zeros(runs+1,inits); beta = zeros(runs+1,inits);
sigma = zeros(2,2,runs+1,inits);
for i = 1:inits
    init(i,:) = [1,1]; initsig(:,:,i)= [1,1;1,2];
    [A,B,S] = gibbs(Y,[X,P],runs,init(i,:),initsig(:,:,i));
    scale = squeeze(S(1,1,:));
end
    beta = B./sqrt(scale);
    alpha = A./sqrt(scale);
    sigma = S./repmat(S(1,1,:),2);

%% Graphs
