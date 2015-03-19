function [] = sim(F, M)
    % SIM - Simulate market equilibrium of entry model.
    %
    %
    % Input arguments:
    %   M = number of markets to simulate
    %   F = number of potential entrants to simulate
    % Outputs:
    %
    %
    
    % model parameters
    MU = 2;
    SIGMA = 1;
    ALPHA = 1;
    BETA = 1;
    DELTA = 1;

    for m = 1:M  % simulate each market individually
        
        Xm = normrnd(3, 1);        % market characteristics
        Zf = normrnd(0, 1, F, 1);  % firm-market characteristics
        
        % fixed costs for each potential entrant
        phi = ALPHA * Zf + normrnd(MU, SIGMA, F, 1);

        profit = @(n) BETA * Xm - DELTA * log(n) - phi;
    
    end

end
