function [mrkts, costs, firms, entry] = sim_markets(F, M)
    % SIM - Simulate market equilibrium of entry model.
    %   Simulate an entry game in M markets with F potential entrants in which
    %   firms make entry decisions sequentially, ordered by fixed costs.
    % Input arguments:
    %   M = number of markets to simulate
    %   F = number of potential entrants to simulate
    % Outputs:
    %   mrkts = M by F+2 matrix with columns for market characteristic,
    %           number of entrants, and profits of each firm
    %   costs = M by F matrix with fixed cost of each firm in columns
    %   firms = M by F matrix with firm-market characteristics in columns
    %   entry = M by F matrix with dummies for entry decisions
    
    % model parameters
    MU = 2; SIGMA = 1;
    ALPHA = 1; BETA = 1; DELTA = 1;

    % allocate space for outcome in each market (simulated below)
    mrkts = zeros(M, F + 2);  % market characteristic, #entrants, and profits
    costs = zeros(M, F);      % fixed costs
    firms = zeros(M, F);      % firm-market characteristics
    entry = zeros(M, F);      % dummies for entry decisions

    for m = 1:M  % simulate each market individually
        
        Xm = normrnd(3, 1);        % market characteristics
        Zf = normrnd(0, 1, F, 1);  % firm-market characteristics
        
        % fixed costs for each potential entrant
        phi = ALPHA * Zf + normrnd(MU, SIGMA, F, 1);
        
        % order firm variables by fixed costs
        [~, idx] = sortrows(phi);
        phi = phi(idx);
        Zf = Zf(idx);

        % profits that would result from sequential entry of each firm
        profits = @(n) BETA * Xm - DELTA * log(n) - phi;
        sequential_profits = profits([1:3]');

        % which firm is the last firm able to earn a non-negative profit?
        entrants = sum(sequential_profits >= 0);
        % what are profits for all firms when this many firms enter?
        realized_profits = profits(entrants);
        realized_profits(entrants + 1:end) = 0;  % non-entrants earn zero
    
        % save this market data for estimation problems
        mrkts(m,:) = [Xm, entrants, realized_profits'];
        costs(m,:) = phi';
        firms(m,:) = Zf';
        entry(m,:) = [ones(1, entrants), zeros(1, F - entrants)];
    end
end
