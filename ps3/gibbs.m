%% Question 2 - Gibbs
% % a) Simulate the dataset

% N=100; K=2; SIGMA=[1,1;1,2]; BETA=1; ALPHA=1;
% [Y, X, P] = sim_dataset(N, k, SIGMA, BETA, ALPHA);
function [ALPHA, BETA, SIGMA] = gibbs(Y,X,runs)
    % Returns results from Gibbs Multinomial Probit
    % INPUTS
    %   Y = N by 1 vector
    %   X = N*K by P vector
    %        P denotes the number of parameters 
    %   runs = number of runs
    % OUTPUTS
    %   ALPHA = runs by 1 vector
    %   BETA = runs by 1 vector
    %   SIGMA = 2*runs by 2 matrix
    global K; K = size(X,2); % K is the number of choice options
    global N; N = 100; 
    
    % Declare outputs
    ALPHA=zeros(runs,1); BETA=zeros(runs,1);
    SIGMA=zeros(2,2,runs); 
    
    % Initialize parameters
    ALPHA(1)=normrnd(0,1); BETA(1)=normrnd(0,1);
    SIGMA(:,:,1)=zeros(2); sigmaold = eye(K);
    SIGMA(:,:,1)=sigmaold;
    if K==2
        betaold = [BETA(1); ALPHA(1)];
    else
        betaold = BETA(1);
    end
    
    % Set priors (Based on R Code by Rossi)
    betabar = zeros(K,1); % Prior mean for BETA
    A = 0.01 * eye(K); % Prior variance for BETA
    NU = K+2; % Prior DoF for SIGMA
    V = NU * eye(K); % Prior scale matrix for SIGMA
    y_init = zeros(K,1);
    
    % Run Gibbs procedure
    for i = 1:runs
        % Draw y_star from truncated normal given ALPHA, BETA, and SIGMA
        [y_star,y_term] = drawystar(X,Y,betaold,sigmaold,y_init);
        
        % Draw ALPHA, BETA given y_star and SIGMA
        betanew = drawbeta(y_star, X, betabar, sigmaold,A);
        %alphanew = drawbeta(y_star, P, betabar, SIGMA,A);
        
        % Draw SIGMA given y_star, ALPHA, BETA
        sigmanew = drawsigma(X, betanew, y_star, V, NU);
        
        % Save results
        ALPHA(i+1) = betanew(2);
        BETA(i+1) = betanew(1);
        SIGMA(:,:,i+1) = sigmanew;
        
        % Update 
        betaold=betanew; sigmaold=sigmanew;
        y_init = y_term;
        
        % Print every 100 iterations
    end
end

function [mu, tau] = postmean(X_i,beta,sigma,w_i,j)
    % Returns posterior mean (2by1) and variance (2by2) from the given data
    % INPUTS
    %   X_i = 2 by 2 matrix of data ([x_i, p_i] in pset)
    %   beta = 2 by 1 vector of coefficients ([BETA, ALPHA] in pset)
    %   sigma = 2 by 2 covariance matrix
    %   w_i = 2 by 1 vector
    %   j = Either 1 or 2.
    % OUTPUTS
    %   mu = p by 1 vector of means
    %   sig = p by p matrix of covariances
    global K;
%     % following book
%     siginv = sigma\eye(K);
%     gamma = siginv(j,:);
%     gamma(j) = [];
%     F = -siginv(j,j) * gamma;
%     xij = X_i(j,:)';
%     xi_j = X_i;
%     xi_j(j,:)=[]; % delete jth row
%     wi_j = w_i;
%     wi_j(j) = [];
%     
%     % Check for empty matrices
%     if isempty(F)
%         F=0;
%     end
%     if isempty(xi_j)
%         xi_j=0;
%     end
%     if isempty(wi_j)
%         wi_j=0;
%     end
%     mu = xij'*beta + F'*(wi_j-xi_j*beta);
%     tau = 1/siginv(j,j);
    % Simply drawing from multivariate normal
    mu = X_i * beta;
    tau = sigma\eye(K);
end

function [y_star, y_term] = drawystar(X, Y, beta, sigma,y_init)
    % Draws y_star from a series of truncated Normal distributions
    % INPUT
    %   X = N*K by p vector of observed data
    %   Y = N by 1 vector of outcomes
    %   beta = p by 1 vector of coefficients
    % OUTPUT
    %   y_star = N*K by 1 vector of imputed latent values
    % Initial value for ystar is zero
    global K; global N;
    N = length(Y);
    
    % Declare variables 
    y_star = zeros(N*K,1); % y_star stores N values of y_star_i
    ystar_old = y_init; % ystar_old is old value of y_star_i; to be replaced by y_star_new
    
    % For each observation, draw over the posterior
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
        %  i=1;j=1  |    i=1;j=2    |    i=2;j=1
        %-------------------------------------------
        %  beta_10  |-> beta_10_new \|  beta_10_new 
        %  beta_11 /|    beta_11     |-> beta_11_new 
        for j=1:K-1
            % Draw from truncated normal
            cont = true;
            c_count = 1;
            while cont
                ystar_new = ystar_old;
                % Calculate posterior conditional means
                z_i = X((i-1)*K+1:i*K,:); % [x_i, p_i]=2by2 matrix
                [mu, tau] = postmean(z_i, beta, sigma,ystar_new,j);
                
                % Draw from posterior using simple acceptance algorithm
%                 ystar_new(j) = normrnd(mu,tau);
                ystar_new = mvnrnd(mu,tau)';
                Ay = weight * ystar_new;
                if y == 1
                    % Need Ay>0
                    cont = not(all(Ay >= 0));
                else
                    % Need Ay <0
                    cont = not(all(Ay <= 0));
                end
                % Update ystar
                
                c_count = c_count+1;
                if mod(c_count,10000) == 0
                    disp(['c_count, y, eval'])
                    disp([c_count,y])
                    disp(weight * ystar_new)
                end
            end
            ystar_old = ystar_new;
        end
        y_star((i-1)*K+1:i*K,:) = ystar_new;
    end
    y_term = ystar_new;
end

function beta = drawbeta(y_star, X, betabar, sigma, A)
    % Draws beta from posterior Normal distribution 
    % INPUTS
    %   y_star = Latent variables that 
    %   SIGMA
    % OUTPUTS
    %   beta =
    % Transform variables to and stack
    global K; global N;
    G = sigma\eye(K); % G is sigma inverted
    C = chol(G,'lower');
    Xstar = kron(eye(N),C')*X;
    ystar = kron(eye(N),C')*y_star;
    % Calculate mean and variance (p.213 McCulloch & Rossi 1994)
    sig = (Xstar' * Xstar+A)\eye(K);
    betahat = sig * (Xstar'*ystar+A*betabar);
    
    % Draw random normal
    beta = mvnrnd(betahat,sig)';
end

function sigma = drawsigma(X,b,ystar,V,NU)
    % Draw sigma from inverted Wishart distribution
    % INPUTS
    %   X = N*K by 2 matrix
    %   b = 2 by 1 vector (=[BETA,ALPHA])
    %   ystar = N*K by 1 vector of latent utilities
    global N; global K;
    % Run regression, get epsilons
    epsilon = reshape(ystar - X*b,2,[]);
%     S = epsilon*epsilon';
    S = diag(sum(epsilon.^2,2));
    
    % Draw from posterior
    % VERSION 1
%     G = wishrnd(V + epsilon'*epsilon,NU+N);
%     sigma = G\eye(K);
    % VERSION 2
    siginv = wishrnd((V + S)\eye(K),NU+N);
    sigma = siginv \ eye(K);
    % VERSION 3 - book version
%     sigma = iwishrnd(V + epsilon*epsilon',NU+N);
end