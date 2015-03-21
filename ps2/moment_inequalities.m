% moment_inequalities - Simulation estimator for entry without
% assumptions on order of entry. Procedure based on Ciliberto and Tamer
% (2009).
    
load data/entry.mat
    
% model parameters (to be converted into function form)
MU = 2; SIGMA = 1;
ALPHA = 1; BETA = 1; DELTA = 1;
NumSims = 100;

% calculate how many firms and markets there are based on inputs
NumMarkets = size(firms, 1);
NumFirms = size(firms,2);
NumConfigs = 2^NumFirms;

% draw u
for i=1:NumSims
    u(:,:,i) = normrnd(MU, SIGMA, NumMarkets, NumFirms);
end

% Get the 2^NumFirms possible market configurations for each market
config = [];
for i = 0:NumFirms
    config = [config; perms([ones(1,i) zeros(1,NumFirms-i)])];
end
config = unique(config, 'rows');

% Loop over all markets and compute profits, entry, and H1 H2
H1hat = [];
H2hat = [];
for m = 1:NumMarkets
% Loop over each simulation and compute profits, entry, and H1 H2
H1hati = [];
H2hati = [];
for i = 1:NumSims
    % Compute profits
    phi = ALPHA*firms(m,:) + u(m,:,i);
    profits_config = repmat(markets(m,1)*BETA - DELTA*log(sum(config,2)),1,NumFirms)- repmat(phi,NumConfigs,1);
    profits_config_plus1 = repmat(markets(m,1)*BETA - DELTA*log(sum(config,2)+1),1,NumFirms)- repmat(phi,NumConfigs,1);

    % Compute entry under each configuration
    entry_config = zeros(NumConfigs, NumFirms);
    for j = 1:NumFirms
        for k = 1:NumConfigs
            if (profits_config(k,j)>=0 & config(k,j)==1) | (profits_config_plus1(k,j)>=0 & config(k,j)==0)
                entry_config(k,j) = 1;
            end
        end
    end
    
    % Check if that entry configuration is an equilibrium
    equi_config = zeros(NumConfigs,1);
    for k = 1:NumConfigs
        if entry_config(k,:)==config(k,:)
            equi_config(k) = 1;
        end
    end
    
    % Compute H1 and H2
    sum_equi_config = sum(equi_config);
    H1 = zeros(NumConfigs,1);
    H2 = zeros(NumConfigs, 1);
    for j = 1:NumConfigs
        if sum_equi_config==1
            H1=equi_config;
            H2=equi_config;
        end
        if sum_equi_config>1
            H2=equi_config;
        end
    end
    H1hati = [H1hati, H1];     
    H2hati = [H2hati, H2];
end
H1hati = mean(H1hati,2);
H2hati = mean(H2hati,2);
H1hat = [H1hat, H1hati];
H2hat = [H2hat, H2hati];
end

H1hat
H2hat