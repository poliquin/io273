function [theta31hat,theta31hatse,theta32hat,theta32hatse] = transitionprob(xt,it)

% Setup some empty matrices
plus1 = zeros(size(xt,1)-1,size(xt,2));
plus2 = zeros(size(xt,1)-1,size(xt,2));
noinvest = zeros(size(xt,1)-1,size(xt,2));
% For every xt where no investment was made, count how many times the state
% went up by 1 and 2. I don't go further because the data did not have the
% state changing by any other number. I also don't investigate the
% investment case since by definition the state following replacement goes
% to 0.
for i=1:(size(xt,1)-1)
    for j=1:size(xt,2)
        if it(i,j)==0
            noinvest(i,j)=1;
            if xt(i,j)+1==xt(i+1,j)
                plus1(i,j)=1;
            elseif xt(i,j)+2==xt(i+1,j)    
                plus2(i,j)=1;
            end
        end
        
    end
end
% Define transition probabilities based on that, also compute standard
% errors for
theta31hat = sum(sum(plus1))/sum(sum(noinvest));
theta31hatse = sqrt(theta31hat*(1-theta31hat)/sum(sum(noinvest)));
theta32hat = sum(sum(plus2))/sum(sum(noinvest));
theta32hatse = sqrt(theta32hat*(1-theta32hat)/sum(sum(noinvest)));
