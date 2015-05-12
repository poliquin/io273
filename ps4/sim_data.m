function [data, it] = sim_data(EV,beta, RC, theta1, theta30, theta31)

data = zeros(1000, 100); % data with 100 columns (buses) and 1,000 rows (time)
it = zeros(1000,100);
for i=1:size(data,2)
    epsilon(:,:,i) = -evrnd(0.5772156649,1,1000,2); % Note: the mu parameter being
    % Euler's constant give a mean of 0. The scale parameter of 1 give a
    % variance of pi^2/2, as required. This can be checked using [M,V] =
    % evstat(0.5772156649,1). Negative is needed because we're modelling
    % maxima.
end
for i=1:size(data,2)
    for j=1:size(data,1)
        [~,choice] = max([-theta1*data(j,i)+epsilon(j,1,i)+beta*EV(data(j,i)+1,1), ...
            -RC-theta1*0+epsilon(j,2,i)+beta*EV(data(j,i)+1,2)]);
        if j<1000 % avoid the boundary error
        if choice==1
            r = rand;
            if r<theta30
                data(j+1,i)=data(j,i);
            elseif r<theta30 + theta31
                data(j+1,i) =data(j,i)+1;
            else
                data(j+1,i) = data(j,i)+2;
            end    
        else
            data(j+1,i)=0;
            it(j,i)=1;
        end
        end
    end
end