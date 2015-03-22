% moment_inequalities - Simulation estimator for entry without
% assumptions on order of entry. Procedure based on Ciliberto and Tamer
% (2009).
    
function obj = moment_inequalities(theta, mu, markets,firms, u)    
% model parameters 
MU = mu; 
ALPHA = theta(1); BETA = theta(2); DELTA = theta(3); SIGMA = theta(4);

% calculate how many firms and markets there are based on inputs
NumMarkets = size(markets, 1);
NumFirms = size(markets,2)-2;
NumConfigs = 2^NumFirms;
NumSims = size(u,3);

% Get the 2^NumFirms possible market configurations for each market
config = [];
for i = 0:NumFirms
    config = [config; perms([ones(1,i) zeros(1,NumFirms-i)])];
end
config = unique(config, 'rows');

% Loop over each market and compute profits, entry, and H1 H2
H1hat = [];
H2hat = [];
for m = 1:NumMarkets
% Loop over each simulation and compute profits, entry, and H1 H2
H1hati = [];
H2hati = [];
for i = 1:NumSims
    % Compute profits
    phi = ALPHA*firms(m,:) + mu + u(m,:,i);
    profits = @(n) repmat(markets(m, 1)*BETA - DELTA*log(n), 1, NumFirms) ...
                   - repmat(phi, NumConfigs, 1);

    profits_config = profits(sum(config, 2));
    profits_config_plus1 = profits(sum(config, 2) + 1);

    % Compute entry under each configuration. Entry occurs if the config
    % calls for the firm to enter and the firm earns positive profits when
    % doing so. It also occurs if the config does not call for the firm to
    % enter but the firm can still earn positive profits by entering after
    % all the firms in the config already entered. If that happens, the
    % config is not a market equilibrium.
    entry_config = zeros(NumConfigs, NumFirms);
    for j = 1:NumFirms
        for k = 1:NumConfigs
            if (profits_config(k, j) >= 0 & config(k, j) == 1) ...
                    | (profits_config_plus1(k, j) >= 0 & config(k, j) == 0)
                entry_config(k, j) = 1;
            end
        end
    end
    
    % Check if that entry configuration is an equilibrium
    equi_config = zeros(NumConfigs,1);
    for k = 1:NumConfigs
        if entry_config(k,:) == config(k,:)
            equi_config(k) = 1;
        end
    end
    
    % Compute H1 and H2
    sum_equi_config = sum(equi_config);
    H1 = zeros(NumConfigs, 1);
    H2 = zeros(NumConfigs, 1);
    for j = 1:NumConfigs
        if sum_equi_config == 1
            H1 = equi_config;
            H2 = equi_config;
        end
        if sum_equi_config > 1
            H2 = equi_config;
        end
    end
    H1hati = [H1hati, H1];     
    H2hati = [H2hati, H2];
end  % end loop over NumSims
H1hati = mean(H1hati, 2);
H2hati = mean(H2hati, 2);
H1hat = [H1hat, H1hati];
H2hat = [H2hat, H2hati];
end  % end loop over NumMarkets

% Determine actual entry
actual_entry = zeros(NumConfigs, NumMarkets);
for mkt = 1:NumMarkets
    entry_temp = markets(mkt, 3:NumFirms+2);
    for firm = 1:NumFirms
        if entry_temp(1,firm) > 0 
            entry_temp(1,firm) = 1;
        end
    end
    for j = 1:NumConfigs
        if entry_temp == config(j,:)
            actual_entry(j, mkt) = 1;
        end
    end
end

% Calculate objective function
% THIS OBJECTIVE FUNCTION IS NOT YET CORRECT - NEED TO REPLACE ACTUAL ENTRY
% WTIH ENTRY PROBABILITIES ESTIMATED NON-PARAMETRICALLY OR A MLOGIT IN STATA.
% AS IS, THOUGH, IT WORKS, BUT IT IS LIKELY NOT GUARANTEED TO BE CONSISTENT.
obj_h1 = actual_entry - H1hat;
obj_h1(obj_h1>0)=0; 
obj_h2 = actual_entry - H2hat;
obj_h2(obj_h2<0)=0; 
obj=[];
for mkt = 1:NumMarkets
    obj = [obj, norm(obj_h1(:,mkt)) + norm(obj_h2(:,mkt))];
end
obj = mean(obj);

