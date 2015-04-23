clear all

rng(8675309);
dir = mfilename('fullpath');
cd(dir(1:end-4));

%% Simulate dataset
N=100;  K=2; SIGMA=[1,1;1,2]; ALPHA=1; BETA=1;
[Y, X, P,As] = sim_dataset(N, K, SIGMA, ALPHA, BETA);

%% Run gibbs using three different starting points and 10000 runs
runs = 10000; inits = 1;
init = zeros(inits,2); initsig = zeros(2,2,inits);
alpha = zeros(runs+1,inits); beta = zeros(runs+1,inits);
sigma = zeros(2,2,runs+1,inits);
for i = 1:inits
    init(i,:) = mvnrnd([1,1],100*eye(2)); initsig(:,:,i)= wishrnd(eye(2),30);
    [A,B,S] = gibbs(Y,[X,P],runs,init(i,:),initsig(:,:,i));
    scale = squeeze(S(1,1,:));
    beta(:,i) = B./scale;
    alpha(:,i) = A./scale;
    sigma(:,:,:,i) = S./repmat(S(1,1,:),2);
end

%% Graphs
