function [prob] = ccp(xt, it)
    % CCP  Estimate conditional choice probabilities.
    %   Find probability of renewal for each mileage state, apply some binning.
    % Input arguments:
    %   xt = matrix of states (mileage of bus in rows, buses in columns)
    %   it = renewal decisions (0 is continue, 1 is renew, buses in columns)
    % Output arguments:
    %   prob = vector of renewal probabilities, indexed by mileage

    states = xt(:);
    choice = it(:);
    counts = tabulate(states);
    N = max(states) + 4;
    
    prob = zeros(N, 1);

    for s=0:N-1
        % get number of times we see this state
        tot = counts(counts(:, 1) == s, 2);
        if isempty(tot)  % we never saw the state
            tot = 0;
            prob(s + 1) = 0;
            continue
        end

        % get number of times we renew in this state
        renewals = sum(choice(states == s));

        % use frequency estimator of probability
        prob(s + 1) = renewals / tot;
    end

    % bin mileages greater than 28
    tot = sum(counts(counts(:, 1) > 28, 2));
    prob(29:end) = sum(choice(states > 28)) / tot;

    % bin mileages less than or equal to 2
    tot = sum(counts(counts(:, 1) <= 2, 2));
    prob(1:2) = sum(choice(states <= 2)) / tot;
end
