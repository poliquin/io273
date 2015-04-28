%% Question 2 - Gibbs
% These functions estimate the Gibbs sampler for pset 3.
function [ALPHA, BETA, SIGMA, ystr] = gibbs(Y,X,runs,initbeta,initsigma)
    % Returns results from Gibbs Multinomial Probit
    % INPUTS
    %   Y = N by 1 vector
    %   X = N*K by P vector
    %        P denotes the number of parameters 
    %   runs = number of runs
    %   initbeta = 1 by K row vector of initial values for [beta, alpha]
    %   initsigma = initial value for sigma 
    % OUTPUTS
    %   ALPHA = runs by 1 vector
    %   BETA = runs by 1 vector
    %   SIGMA = 2*runs by 2 matrix
    global sigsize; sigsize = length(initsigma);
    global K; K=2;% K is the number of choice options
    global N; N = 100; 
    
    % Allocate space for outputs
    ALPHA=zeros(runs,1); BETA=zeros(runs,1);
    SIGMA=zeros(sigsize,sigsize,runs); 
    ystr=zeros(N,sigsize,runs);
    
    % Initialize parameters
    ALPHA(1)=initbeta(2); BETA(1)=initbeta(1);
    SIGMA(:,:,1)=initsigma; sigmaold = SIGMA(:,:,1);
    betaold = [BETA(1); ALPHA(1)];
   
    % Set diffuse priors (Based on R Code by Rossi)
    betabar = zeros(K,1); % Prior mean for BETA
    A = 0.01 * eye(K); % Prior inverse variance for BETA
    NU = sigsize+2; % Prior DoF for SIGMA
    V = NU * eye(sigsize); % Prior scale matrix for SIGMA
    
    P=reshape(Y(:,2:end)',[],1);
    Z=X(:,2);
    
    % Run Gibbs procedure
    for i = 1:runs
        % 1) Draw y_star (Nx2) from a posterior Normal
        % 2) Draw beta (2x1) from a posterior Normal
        % 3) Draw sigma from an inverse Wishart
%         if mod(i,100)==0
%             fprintf('i %i\n',i)
%             fprintf('beta %3.3f, gamma %3.3f\n',betaold(1),betaold(2))
%             disp(sigmaold)
%         end
        [y_star] = drawystar(X,Y,betaold,sigmaold);
%         ystr(:,:,i)=y_star;
        betanew = drawbeta(y_star, X, betabar, sigmaold,A);
        sigmanew = drawsigma(Y,X, betanew, y_star, V, NU);
        
        % Save results
        ALPHA(i+1) = betanew(2);
        BETA(i+1) = betanew(1);
        SIGMA(:,:,i+1) = sigmanew;
        
        % Update y*, beta, sigma
        betaold=betanew; sigmaold=sigmanew;
        
        % Print every 1000 iterations
        if mod(i,1000)==0
            fprintf('Iteration i = %i\n',i)
        end
    end
end

function [y_star] = drawystar(X, Y, beta, sigma)
    % Draws y_star from a Normal distribution centered at Xbeta with
    % variance sigma
    % INPUT
    %   X = N*K by p vector of observed data
    %   Y = N by K vector of outcomes
    %   beta = p by 1 vector of coefficients. Always beta(1) = beta with 
    %           beta(2) = alpha or gamma.
    %   sigma = covariane matrix
    %   y_init = initial value for y_star
    % OUTPUT
    %   y_star = N by K vector of imputed latent values (y1*, y2*)
    %   y_term = ending value of y_star; use for next iteration.
    
    % Initial value for ystar is zero
    global sigsize; global N; global K;
    
    % Declare variables 
    y_star = zeros(N*2,1); % y_star stores N values of y_star_i
    if sigsize == 4
        P=reshape(Y(:,2:end)',[],1);
        % Need to draw eta
        eta = P-beta(2)*X(:,2);
%         tempgam = (X(:,2)'*X(:,2))\(X(:,2)'*P);
%         eta= P - tempgam*X(:,2); % eta = P-Zgamma (2N by 1)
        sig11 = sigma(1:2,1:2); sig22 = sigma(3:4,3:4);
        sig12 = sigma(1:2,3:4); sig21 = sigma(3:4,1:2);
        sigma = sig11-sig12*(sig22\sig21);
    end
    
    % For each observation, draw over the posterior
    for ind=1:N
        % Initialize some parameters
        y = Y(ind);
        weight = cmat(ind);
        
        % Update beta_i
        cont = true;
        c_count = 1;
        while cont
            % Calculate posterior conditional means
            x_i = X((ind-1)*K+1:ind*K,:);
            if sigsize==4
                % need to draw from y* ~ N(Xbeta + P + condmean(epsilon),
                % condvar(epsilon))
                mu = x_i(:,1)*beta(1) + P((ind-1)*K+1:ind*K,:)+ ... 
                    sig12*(sig22\eta((ind-1)*K+1:ind*K));
                tau = sigma;
%                 T = [1,0,1,0;0,1,0,1];
%                 mu = x_i*beta;
%                 tau = T*sigma*T';
            else
                mu=x_i*beta;
                tau = sigma;
            end
            
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

            c_count = c_count+1;
            if mod(c_count,100000) == 0
                % If it gets stuck, push it a little bit
                sigma = sigma * 100;
                disp('Got stuck; multiply sigma by 100')
                disp(['i, c_count, y, eval'])
                disp([ind,c_count,y])
                disp(weight * ystar_new)
                disp(sigma)
            end
        end
        % Update ystar
        y_star((ind-1)*K+1:ind*K,:) = ystar_new;
    end
end

function beta = drawbeta(y_star, X, betabar, sigma, A)
    % Draws beta from posterior Normal distribution centered at 
    % INPUTS
    %   y_star = N*K by 1 matrix of Latent variables
    %   sigma = K by K covariance matrix
    %   X = [X, P] in pset 3
    %   betabar = prior mean on beta
    %   sigma = prior on sigma
    %   A = prior variance on beta
    % OUTPUTS
    %   beta = K by 1 posterior draw of beta
    % Transform variables by multipyling by C and stack
    global K; global N; global sigsize;
    if sigsize ==4
        % Change variance matrix to correct dimensions.
        % Is this right??
        T=[1,0,1,0;0,1,0,1];
        sigma = T*sigma*T';
%         % Or is this right??
%         sig11 = sigma(1:2,1:2); sig22=sigma(3:4,3:4);
%         sig12 = sigma(1:2,3:4); sig21 = sigma(3:4,1:2);
%         sigma = sig11-sig12*(sig22\sig21);
    end
    G = sigma\eye(2); % G is sigma inverted
    C = chol(G,'lower');
    Xstar = kron(eye(N),C)*X;
    ystar = kron(eye(N),C)*y_star;
    sig = (Xstar'*Xstar + A)\eye(K);
    betahat = sig * (Xstar'*ystar+A*betabar);
    
    % Draw random normal
    beta = mvnrnd(betahat,sig)';
end

function sigma = drawsigma(Y,X,b,ystar,V,NU)
    % Draw sigma from inverted Wishart distribution
    % INPUTS
    %   Y = [Y, P1, P2] matrix (N by 3)
    %   X = N*K by 2 matrix of either [X,P] or [X,Z]
    %   b = 2 by 1 vector (=[BETA,ALPHA])
    %   ystar = N by K vector of latent utilities
    %   NU = scalar addition to dof in IW distribution
    % OUTPUTS
    %   sigma = sigsize * sigsize covariance matrix
    global N; global sigsize;
    % Run regression, get epsilons
    if sigsize == 2
        error=reshape(ystar-X*b,2,[])'; %N by 2
    elseif sigsize == 4
        Z = X(:,2);
        X = X(:,1);
        P = reshape(Y(:,2:3)',[],1);
        beta = b(1);
        gamma = b(2);
        etahat=P-Z*gamma; % 2N by 1
        epsilon=ystar-X*beta-P; % 2N by 1
%         etahat = draweta(P, Z, gamma, epsilon, omega);
        % Reshape
        eta = reshape(etahat,2,[])';
        epsilon=reshape(epsilon,2,[])';
        error = [epsilon,eta]; % N by 4
    end
    % Draw from posterior
    sigma = iwishrnd(V + error'*error,NU+N);
end
% 
% function etahat = draweta(P, Z, gam, eps, sigma)
%     % Draw eta by Gibbs'ing through the observations from posterior 
%     % conditional posterior normal centered at P-Z*gamma, given epsilon
%     % INPUTS
%     %   P = N*K vector 
%     %   Z = N*K vector
%     %   gam = scalar parameter
%     %   epsilon = N*K vector of imputed errors
%     %   Omega = sigsize*sigsize matrix of covariances
%     % OUTPUTS
%     %   etahat = N*K vector of draw from conditional posterior distribution
%     global N; global K;
%     etahat = zeros(N*K,1);
%     % Get covariance matrix and mean conditional on epsilon
%     sig11 = sigma(1:2,1:2); sig22=sigma(3:4,3:4);
%     sig12 = sigma(1:2,3:4); sig21 = sigma(3:4,1:2);
%     postsig = sig22 - sig21 * (sig11\sig12);
%     % This is the mean? Maybe this should be zero
%     tempgam = (Z'*Z)\(Z'*P);
%     eta = P - Z*tempgam;
% %     etahat_i = zeros(K,1); % Initial value for etahat
%     for ind = 1:N
%         % Init variables 
%         epsilon_i = eps((ind-1)*K+1:ind*K);
%         etaold = eta((ind-1)*K+1:ind*K);
% %         etaold = etahat_i; %update old value
%         % Get posterior mean centered on old value
% %         postmu = etaold + sig21 * (sig11\epsilon_i);
%         postmu = sig21 * (sig11\epsilon_i);
%         % Draw from posterior and save
%         etahat_i = mvnrnd(postmu,postsig)';
%         etahat((ind-1)*K+1:ind*K)=etahat_i;
%     end
% end

function weight = cmat(i)
    if mod(i,2)==0
        weight = [-1,1;1,1];
    else
        weight = [1,1;1,-1];
    end
end
