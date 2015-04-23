%% Question 2 - Gibbs
% These functions estimate the Gibbs sampler for pset 3.
function [ALPHA, BETA, SIGMA] = gibbs(Y,X,runs,initbeta,initsigma)
    % Returns results from Gibbs Multinomial Probit
    % INPUTS
    %   Y = N by 1 vector
    %   X = N*K by P vector
    %        P denotes the number of parameters 
    %   runs = number of runs
    %   initbeta = 1 by K row vector of initial values for [beta, alpha]
    %   initsigma = initial value for sigma (scalar or K by K)
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
    ALPHA(1)=initbeta(2); BETA(1)=initbeta(1);
    SIGMA(:,:,1)=eye(K)*initsigma; sigmaold = SIGMA(:,:,1);
    if K==2
        betaold = [BETA(1); ALPHA(1)];
    else
        betaold = BETA(1);
    end
    
    % Set diffuse priors (Based on R Code by Rossi)
    betabar = zeros(K,1); % Prior mean for BETA
    A = 0.01 * eye(K); % Prior variance for BETA
%     betabar = [2;2];
%     A=eye(K);
    NU = K+2; % Prior DoF for SIGMA
    V = NU * eye(K); % Prior scale matrix for SIGMA
    y_init = zeros(K,1);
    
    % Run Gibbs procedure
    for i = 1:runs
        % 1) Draw y_star from a posterior Normal
        % 2) Draw beta from a posterior Normal
        % 3) Draw sigma from an inverse Wishart
        [y_star,y_term] = drawystar(X,Y,betaold,sigmaold,y_init);
        betanew = drawbeta(y_star, X, betabar, sigmaold,A);
        sigmanew = drawsigma(X, betanew, y_star, V, NU);
        
        % Save results
        ALPHA(i+1) = betanew(2);
        BETA(i+1) = betanew(1);
        SIGMA(:,:,i+1) = sigmanew;
        
        % Update 
        betaold=betanew; sigmaold=sigmanew;
        y_init = y_term;
        
        % Print every 100 iterations
        if mod(i,100)==0
            disp('i = ')
            disp(i)
        end
    end
end

function [y_star, y_term] = drawystar(X, Y, beta, sigma,y_init)
    % Draws y_star from a Normal distribution centered at Xbeta with
    % variance sigma
    % INPUT
    %   X = N*K by p vector of observed data
    %   Y = N by 1 vector of outcomes
    %   beta = p by 1 vector of coefficients
    %   sigma = covariane matrix
    %   y_init = initial value for y_star
    % OUTPUT
    %   y_star = N*K by 1 vector of imputed latent values
    %   y_term = ending value of y_star; use for next iteration.
    
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
        
    % Update beta_i
        % Draw from truncated normal
        cont = true;
        c_count = 1;
        while cont
            % Calculate posterior conditional means
            z_i = X((i-1)*K+1:i*K,:); % [x_i, p_i]=2by2 matrix
            mu = z_i * beta;
            tau = sigma;

            % Draw from posterior using simple acceptance algorithm
            ystar_new = mvnrnd(mu,tau)';
            Ay = weight * ystar_new;
            Ay(Ay<0)=0;
            Ay=all(Ay);
            if y == Ay
                cont = false;
            else
                cont = true;
            end
            % Update ystar

            c_count = c_count+1;
            if mod(c_count,100000) == 0
                % If it gets stuck, push it a little bit
                sigma = sigma * 100;
                disp('Got stuck; multiply sigma by 100')
                disp(['i, c_count, y, eval'])
                disp([i,c_count,y])
                disp(weight * ystar_new)
            end
        end
        ystar_old = ystar_new;
        y_star((i-1)*K+1:i*K,:) = ystar_new;
    end
    y_term = ystar_new;
end

function beta = drawbeta(y_star, X, betabar, sigma, A)
    % Draws beta from posterior Normal distribution centered at 
    % INPUTS
    %   y_star = Latent variables
    %   sigma = K by K covariance matrix
    %   X = [X, P] in pset 3
    %   betabar = prior mean on beta
    %   sigma = prior on sigma
    %   A = prior variance on beta
    % OUTPUTS
    %   beta = K by 1 posterior draw of beta
    % Transform variables by multipyling by C and stack
    global K; global N;
    G = sigma\eye(K); % G is sigma inverted
    C = chol(G,'lower');
    Xstar = kron(eye(N),C)*X;
    ystar = kron(eye(N),C)*y_star;
    % Calculate mean and variance (p.213 McCulloch & Rossi 1994)
    sig = (Xstar'*Xstar + A)\eye(K);
    betahat = sig * (Xstar'*ystar + A*betabar);
    
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
    
    % Draw from posterior
    sigma = iwishrnd(V + epsilon*epsilon',NU+N);
end

% function [mu, tau] = postmean(X_i,beta,sigma,w_i,j)
%     % Returns posterior mean (2by1) and variance (2by2) from the given data
%     % INPUTS
%     %   X_i = 2 by 2 matrix of data ([x_i, p_i] in pset)
%     %   beta = 2 by 1 vector of coefficients ([BETA, ALPHA] in pset)
%     %   sigma = 2 by 2 covariance matrix
%     %   w_i = 2 by 1 vector
%     %   j = Either 1 or 2.
%     % OUTPUTS
%     %   mu = p by 1 vector of means
%     %   sig = p by p matrix of covariances
%     global K;
% %     % following book
% %     siginv = sigma\eye(K);
% %     gamma = siginv(j,:);
% %     gamma(j) = [];
% %     F = -siginv(j,j) * gamma;
% %     xij = X_i(j,:)';
% %     xi_j = X_i;
% %     xi_j(j,:)=[]; % delete jth row
% %     wi_j = w_i;
% %     wi_j(j) = [];
% %     
% %     % Check for empty matrices
% %     if isempty(F)
% %         F=0;
% %     end
% %     if isempty(xi_j)
% %         xi_j=0;
% %     end
% %     if isempty(wi_j)
% %         wi_j=0;
% %     end
% %     mu = xij'*beta + F'*(wi_j-xi_j*beta);
% %     tau = 1/siginv(j,j);
%     % Simply drawing from multivariate normal
%     mu = X_i * beta;
%     tau = sigma;
% end