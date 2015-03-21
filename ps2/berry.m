function [like] = berry(mrkts, firms, entry, mu, sigma, theta, draws)
    % BERRY - Simulation estimator for sequential entry model
    %   Produce estimate in equation (A.2) using equations (12) and (A.1) plus
    %   the assumption that firms enter in order of profitability.
    % Input arguments:
    %   mrkts = M by . matrix with market characteristic in first column
    %   firms = M by F matrix with firm-market characteristics in columns
    %   entry = M by F matrix of dummies for entry decisions
    %   mu = mean of distribution for the error term, u_fm
    %   sigma = standard deviation of distribution for the error term, u_fm
    %   draws = T by F*M matrix of normal draws normal to use in simulation
    %   theta = values for ALPHA, BETA, and DELTA
    % Outputs:
    %   like = estimated log likelihood

    ALPHA = theta(1); BETA = theta(2); DELTA = theta(3);
    [M, F] = size(firms);  % number of markets and potential entrants
    Phat = zeros(M, F);    % hold simulated probabilities of entry
    T = size(draws, 1);    % number of simulation draws

    for m = 1:M  % simulate each market individually

        Xm = mrkts(m, 1);           % market characteristic
        Zf = firms(m, :)';          % firm-market characteristics
        predictions = zeros(F, 1);  % probability of entry for each firm in m
        
        scol = (m - 1) * F + 1;   % first simulation column for this market
        ecol = m * F;             % last simulation column for this market
        
        for t = 1:T  % loop over T draws of simulation estimator

            % simulated fixed costs for each firm, with u = mu + sigma * v
            simphi = ALPHA * Zf + mu + sigma * draws(t, scol:ecol)';
            
            % order firms by simulated fixed costs
            [~, idx] = sortrows(simphi);
            
            % estimated profits that result from entry by n firms
            phat = @(n) BETA * Xm - DELTA * log(n) - simphi;
            % expected number of entrants when there are K possible entrants
            nhat = @(K) sum(arrayfun(@(n) sum(phat(n) >= 0) >= n, 1:K));

            % when firms are ordered by profits, predicted entry equals 1 for
            % firms with rank less than or equal to the expected # of entrants
            predictions = predictions + (1/T) * (idx <= nhat(F));

        end  % end simulation loop

        % estimated probability of entry for each firm in market m
        Phat(m, :) = predictions';
    end  % end market loop
    
    % claculate the log likelihood
    mkt_likelihood = prod((Phat.^entry) .* ((1 - Phat).^(1 - entry)), 2);
    like = sum(log(max(10e-20, mkt_likelihood)));    
end
