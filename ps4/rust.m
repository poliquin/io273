function loglike = rust(theta1, theta31, theta32, RC, beta,xt,it)
% compute the 0 transition probability
theta30 = 1 - theta31 - theta32;
% compute EV function (fixed point iteration)
evhat = ev(.99, 10, theta1, theta30, theta31, theta32);
% setup empty matrix
phat = zeros(size(xt,1),size(xt,2));
% loop over to get probability of investment
for i=1:size(xt,1)
    for j=1:size(xt,2)
        x = xt(i,j);
phat(i,j) = exp(-RC-theta1*0+beta*evhat(1,2))/(exp(-RC-theta1*0+beta*evhat(1,2)) + exp(-theta1*x+beta*evhat(x+1,1)));
    end
end
% define log-likelihood
like = prod((phat.^it).*(1-phat).^(1-it),2);
loglike = sum(log(max(10e-20, like)));    
