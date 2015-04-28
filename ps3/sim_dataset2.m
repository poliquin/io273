function [Y, X, P, Z] = sim_dataset2(N, k, BETA, GAMMA)
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
    Z = zeros(N,k);
    for i = 1:N
        % Draw X, Z, errors.
        OMEGA = [1,1,1,1;1,2,2,2;1,2,3,3;1,2,3,4];
        errors = mvnrnd(zeros(k*2,1),OMEGA);
        epsilon = errors(1:2);
        eta = errors(3:4);
        X_i = mvnrnd(zeros(k,1), eye(k));
        Z_i = mvnrnd(zeros(k,1), eye(k));
        p_i = Z_i * GAMMA+ eta;

        % Calculate y_star (1 by k)
        y_star = X_i * BETA + Z_i * GAMMA + eta + epsilon;

        % Calculate y
        if mod(i,2) == 0
            y = all(A_EVEN*y_star' > 0);
            A(:, :, i) = A_EVEN;
        else
            y = all(A_ODD*y_star' > 0);
            A(:, :, i) = A_ODD;
        end
        Y(i,:) = y;
        X(i,:) = X_i';
        P(i,:) = p_i';
        Z(i,:) = Z_i';
    end
end