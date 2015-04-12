function [Y, X] = sim_dataset(N, k, SIGMA, BETA, ALPHA)
    % Generate dataset for part 2 of the problem set
    % Input Arguments:
    %   N = Number of observations
    %   k = length of y vector
    %   SIGMA = 2 by 2 covariance matrix
    %   BETA = scalar
    %   ALPHA = scalar
    % Outputs:
    %   Y = N by k matrix
    %   X = N by 1 vector
    A_ODD = [1, 1; 1, -1]; A_EVEN = [-1, 1; 1, 1];
    Y = zeros(N,1); X = zeros(N,2);
    for i = 1:N
        % Draw X, p, EPSILON (2 by 1)
        Xi = mvnrnd(zeros(2,1), eye(2));
        p = mvnrnd(zeros(2,1), eye(2));
        EPSILON = mvnrnd(zeros(2,1),SIGMA);

        % Calculate y_star (k by 1)
        y_star = Xi * BETA + p * ALPHA + EPSILON;

        % Calculate y
        if mod(i,2) == 0
            y = all(A_EVEN*y_star' > 0);
        else
            y = all(A_ODD*y_star' > 0);
        end
        Y(i,:) = y;
        X(i,:) = Xi';
    end