function pi = fppi(EV, beta, RC, theta1, theta30, theta31, theta32)

%% Define parameters
TOL     = 1e-9; % tolerance 

%% Define grid and initial guess
pi       = zeros(200,2); 
pi(1,1) = 1;

%% Value function iteration

diff = 1;
while diff > TOL % run this loop until diff drops below TOL
    piold = pi; % save old value function
    
    % construct new value function v
    for j = 1:200 % have to start at 1 because Matlab won't take 0...
        for i = 1:2
            if j==1
            if i==1 % no investment
                pi(j,i) = (sum(piold(:,2))+piold(1,1)*theta30)*exp(beta*EV(j,1))/(exp(beta*EV(j,1))+exp(-RC+beta*EV(1,2)));
            else % investment
                pi(j,i) = (sum(piold(:,2))+piold(1,1)*theta30)*exp(-RC+beta*EV(1,2))/(exp(beta*EV(j,1))+exp(-RC+beta*EV(1,2)));                
            end
            elseif j==2
            if i==1 % no investment
                pi(j,i) = piold(j,1)*theta30*exp(-theta1*(j-1)+beta*EV(j,1))/(exp(-theta1*(j-1)+beta*EV(j,1))+exp(-RC+beta*EV(1,2))) ...
                + piold(j-1,1)*theta31*exp(-theta1*(j-2)+beta*EV(j-1,1))/(exp(-theta1*(j-2)+beta*EV(j-1,1))+exp(-RC+beta*EV(1,2)));
            else % investment
                pi(j,i) = piold(j,1)*theta30*exp(-RC+beta*EV(1,2))/(exp(-theta1*(j-1)+beta*EV(j,1))+exp(-RC+beta*EV(1,2))) ...
                + piold(j-1,1)*theta31*exp(-RC+beta*EV(1,2))/(exp(-theta1*(j-2)+beta*EV(j-1,1))+exp(-RC+beta*EV(1,2)));
            end
            elseif j>2
            if i==1 % no investment
                pi(j,i) = piold(j,1)*theta30*exp(-theta1*(j-1)+beta*EV(j,1))/(exp(-theta1*(j-1)+beta*EV(j,1))+exp(-RC+beta*EV(1,2))) ...
                + piold(j-1,1)*theta31*exp(-theta1*(j-2)+beta*EV(j-1,1))/(exp(-theta1*(j-2)+beta*EV(j-1,1))+exp(-RC+beta*EV(1,2))) ...
                + piold(j-2,1)*theta32*exp(-theta1*(j-3)+beta*EV(j-2,1))/(exp(-theta1*(j-3)+beta*EV(j-2,1))+exp(-RC+beta*EV(1,2)));
            else % investment
                pi(j,i) = piold(j,1)*theta30*exp(-RC+beta*EV(1,2))/(exp(-theta1*(j-1)+beta*EV(j,1))+exp(-RC+beta*EV(1,2))) ...
                + piold(j-1,1)*theta31*exp(-RC+beta*EV(1,2))/(exp(-theta1*(j-2)+beta*EV(j-1,1))+exp(-RC+beta*EV(1,2))) ...
                + piold(j-2,1)*theta32*exp(-RC+beta*EV(1,2))/(exp(-theta1*(j-3)+beta*EV(j-2,1))+exp(-RC+beta*EV(1,2)));
            end
            end
        end
    end
    
    diff = max(max(abs(pi-piold)));
end
