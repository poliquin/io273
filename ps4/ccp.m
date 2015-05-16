function [prob] = ccp(xt, it)
    % CCP  Estimate conditional choice probabilities.
    %   Find probability of renewal for each mileage state between 0 and 42.
    % Input arguments:
    %   xt = matrix of states (mileage of bus in rows, buses in columns)
    %   it = renewal decisions (0 is continue, 1 is renew, buses in columns)
    % Output arguments:
    %   prob = vector of renewal probabilities, indexed by mileage

    states = xt(:);
    choice = it(:);
    counts = tabulate(states);

    prob = zeros(max(max(xt))+4, 1);

    for s=0:(max(max(xt))+3)
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

    % bin probs less than or equal to 2
    i=1;
    while prob(1)==0
        i=i+1;
        tot = sum(counts(counts(:, 1) <= i, 2));
        prob(1:i) = sum(choice(states <=i)) / tot;
    end
    
end
