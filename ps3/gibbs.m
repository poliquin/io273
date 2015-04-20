%% Question 2 - Gibbs
% % a) Simulate the dataset

% N=100; K=2; SIGMA=[1,1;1,2]; BETA=1; ALPHA=1;
% [Y, X, P] = sim_dataset(N, k, SIGMA, BETA, ALPHA);
function [ALPHA, BETA, SIGMA] = gibbs(Y,X,runs)
    % Returns results from Gibbs Multinomial Probit
    % INPUTS
    %   Y = N by 1 vector
    %   X = N by K vector
    %   weight = K by K matrix
    %   runs = number of runs
    % OUTPUTS
    %   ALPHA = runs by 1 vector
    %   BETA = runs by 1 vector
    %   SIGMA = 2*runs by 2 matrix
    global K
    K = 2;
    
    % Declare outputs
    ALPHA=zeros(runs,1); BETA=zeros(runs,1);
    SIGMA=zeros(2,2,runs); 
    
    % Initialize parameters
    ALPHA(1)=normrnd(0,1); BETA(1)=normrnd(0,1);
    SIGMA(:,:,1)=zeros(2);
    
    % Set priors (Based on R Code by Rossi)
    betabar = zeros(K,1); % Prior mean for BETA
    A = 0.01 * eye(K); % Prior variance for BETA
    NU = K+2; % Prior DoF for SIGMA
    V = NU * eye(K); % Prior scale matrix for SIGMA
    
    % Run Gibbs procedure
    for i = 1:runs
        % Draw y_star from truncated normal given ALPHA, BETA, and SIGMA
        y_star = drawystar(mu,sig,Y);
        
        % Draw ALPHA, BETA given y_star and SIGMA
        betanew = drawbeta(y_star, X, betabar, SIGMA,A,K,N);
        alphanew = drawbeta(y_star, P, betabar, SIGMA,A,K,N);
        
        % Draw SIGMA given y_star, ALPHA, BETA
        sigmanew = drawsigma();
        
        % Save results
        ALPHA(i+1) = alphanew;
        BETA(i+1) = betanew;
        SIGMA(:,:,i+1) = sigmanew;
        
        % Print every 100 iterations
    end
end

function y_star = drawystar(mu, sig, Y)
    % Draws y_star from a series of truncated Normal distributions
    % INPUT
    %   mu = N by 1 vector of means
    %   sig = N by 1 vector of sd
    %   Y = N by 1 vector of outcomes
    % OUTPUT
    %   y_star = N*K by 1 vector of imputed latent values
    % Initial value for ystar is zero
    N = length(Y);
    y_star = zeros(N*K,1);
    ystar_old = zeros(2,1);
    
    % For each observation
    for i=1:N
        % Initialize some parameters
        y = Y(i);
        % Constraint matrix (matrix A in pset)
        if mod(i,2)==0
            weight = [-1,1;1,1];
        else
            weight = [1,1;1,-1];
        end
        
        % Update beta_i as so
        %  i=1;j=1  |     i=1;j=2   |    i=2;j=1
        %---------------------------------------
        % |beta_10| | |beta_10_new| | |beta_10_new|
        % |beta_11| | |  beta_11  | | |beta_11_new|
        for j=1:2
            % Draw from truncated normal
            cont = true;
            while cont
                % Draw normals until they match data
                ystar_new = ystar_old;
                ystar_new(j) = normrnd(mu,sig);
                % Check whether falls into range
                if y == 1
                    % Need Ay>0
                    cont = weight * ystar_new <= 0;
                else
                    % Need Ay <0
                    cont = weight * ystar_new >= 0;
                end
                % Update ystar
                ystar_old = ystar_new;
            end
        end
        y_star((i-1)*k+1:i*k,:) = ystar_new;
    end
    y_star = ystar_new;
end

function beta = drawbeta(y_star, X, betabar, SIGMA, A, K, N)
    % Draws beta from posterior Normal distribution 
    % INPUTS
    %   y_star = Latent variables that 
    %   SIGMA
    % OUTPUTS
    %   beta =
    % Transform variables
    G = SIGMA\eye(K); % G is sigma inverted
    C = chol(G,'lower');
    % Stack X and Y
    Xstar = kron(eye(N),C')*X;
    ystar = kron(eye(N),C')*y_star;
    % Calculate mean and variance (p.213 McCulloch & Rossi 1994)
    sig = 1/(Xstar' * Xstar+A);
    betahat = sig * (Xstar'*ystar+A*betabar);
    
    % Draw random normal
    beta = normrnd(betahat,sig);
end

function sigma = drawsigma(X,b,ystar)
    % Draw sigma from inverted Wishart distribution
    % INPUTS
    %   X = N*K by 2 matrix
    %   b = 2 by 1 vector (=[BETA,ALPHA])
    %   ystar = N*K by 1 vector of latent utilities
    
    % Run regression, get epsilons
    epsilon = ystar - X*b;
    
    % 
    vtilde = wishrnd(V + 
end