function [v] = ev(beta, RC, theta1, theta30, theta31, theta32)

%% Define parameters
TOL     = 1e-9; % tolerance 

%% Define grid and initial guess
v       = zeros(303,2); 

%% Value function iteration

diff = 1;
while diff > TOL % run this loop until diff drops below TOL
    vold = v; % save old value function
    
    % construct new value function v
    for j = 1:301 % have to start at 1 because Matlab won't take 0...
        for i = 1:2
            if i==1 % no investment
                x=j; 
            else % investment
                x=1;
            end
            v(j,i) = theta30*log(exp(-RC-theta1*0+beta*v(1,2)) + ...
                exp(-theta1*(x-1)+beta*v(x,1))) ...
                + theta31*log(exp(-RC-theta1*0+beta*v(1,2)) + ...
                exp(-theta1*(x)+beta*v(x+1,1)))...
                + theta32*log(exp(-RC-theta1*0+beta*v(1,2)) + ...
                exp(-theta1*(x+1)+beta*v(x+2,1)));
        end
    end
    
    diff = max(abs(v-vold));
end

v = v(1:201,:);
