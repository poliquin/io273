function [Y, X, P, A] = sim_dataset(N, k, SIGMA, BETA, ALPHA)
    % Generate dataset for part 2 of the problem set
    % Input Arguments:
    %   N = Number of observations
    %   k = length of y vector
    %   SIGMA = 2 by 2 covariance matrix
    %   BETA = scalar
    %   ALPHA = scalar
    % Outputs:
    %   Y = N by 1 vector
    %   X = NK by 1 vector
    A_ODD = [1, 1; 1, -1]; A_EVEN = [-1, 1; 1, 1];
    Y = zeros(N,1); X = zeros(N,k); P = zeros(N,k);
    for i = 1:N
        % Draw X, p, EPSILON (2 by 1)
        X_i = mvnrnd(zeros(k,1), eye(k));
        p_i = mvnrnd(zeros(k,1), eye(k));
        EPSILON = mvnrnd(zeros(2,1),SIGMA);

        % Calculate y_star (k by 1)
        y_star = X_i * BETA + p_i * ALPHA + EPSILON;

        % Calculate y
        if mod(i,2) == 0
            y = all(A_EVEN*y_star' > 0);
            A(:, :, i) = A_EVEN;
        else
            y = all(A_ODD*y_star' > 0);
            A(:, :, i) = A_ODD;
        end
        Y(i,:) = y';
        X(i,:) = X_i';
        P(i,:) = p_i';
    end
end