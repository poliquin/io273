function [like] = berry(mrkts, firms, entry, mu, sigma, theta, draws, order)
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
    %   order = 'ascend' means order of entry is most to least profitable,
    %           'descend' means order of entry is least to most profitable
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
            
            % estimated profits that result from entry by n firms
            phat = @(n) BETA * Xm - DELTA * log(n) - simphi;
            % expected number of entrants when there are K possible entrants
            nhat = @(K) sum(arrayfun(@(n) sum(phat(n) >= 0) >= n, 1:K));
            entrants = nhat(F);

            % order firms by simulated fixed costs
            if strcmp(order, 'ascend')  % most profitable enter first
                [~, idx] = sortrows(simphi, 1);
                choices = (idx <= entrants);  % dummy marker for entry
            elseif strcmp(order, 'descend')  % least profitable enter first
                [~, idx] = sortrows(simphi, -1);
                if entrants < F  % need to decide which firms enter
                    % firms that make money in a n + 1 equilibrium surely enter
                    profits_n_plus_one = phat(entrants + 1);
                    choices = profits_n_plus_one >= 0;
                    num_entry = sum(choices);
                    if num_entry < entrants
                        % n-L profitable firms in an n firm equilibrium enter
                        profits_n = phat(entrants);
                        % which firms now wish they could enter?
                        want_enter = profits_n >= 0 & profits_n_plus_one < 0;
                        cs = cumsum(want_enter);
                        % which of the remaining firms get to enter?
                        entry_order = sortrows([idx, cs]);
                        addentry = entry_order(:, 2) <= (entrants - num_entry);
                        new_choices = addentry .* entry_order(:, 2);
                        % create a single vector of entry decisions
                        choices = choices + new_choices(idx);
                    end
                else  % everyone enters, so no need to pick firms
                    choices = ones(3, 1);
                end
            end
            % when firms are ordered by profits, predicted entry equals 1 for
            % firms with rank less than or equal to the expected # of entrants
            predictions = predictions + (1/T) * choices;

        end  % end simulation loop

        % estimated probability of entry for each firm in market m
        Phat(m, :) = predictions';
    end  % end market loop
    
    % claculate the log likelihood
    mkt_likelihood = prod((Phat.^entry) .* ((1 - Phat).^(1 - entry)), 2);
    like = sum(log(max(10e-20, mkt_likelihood)));    
end
