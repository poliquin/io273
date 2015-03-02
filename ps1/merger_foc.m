function foc = merger_foc(P,MC,X,BETA,ALPHA,XI,NU,SIGMA,merger)
    [shares, sim_DS,sim_XDS]=merger_shr(P,X,BETA,ALPHA,XI,NU,SIGMA);
    %MERGER_EQUILIBRIUM: 
    coeffmat = [0,1,0;1,0,1;0,0,0];
    if merger
        foc = shares + (P-MC).*sim_DS + coeffmat*(P-MC).*sim_XDS;
    else 
        foc = shares + (P-MC).*sim_DS;
    end
end

function [shares, sim_DS,sim_XDS]=merger_shr(P0,X,BETA,ALPHA,XI,NU,SIGMA)
    deltas = X * BETA - ALPHA * P0+ XI;
    price_utility = bsxfun(@times,P0,NU);
    f_sim = deltashares(deltas, price_utility, 3);
    shares = reshape(mean(f_sim,2),3,[]);

    % Compute cross prive derivatives, reshape
    sim_DS = -mean(bsxfun(@times,f_sim .* (1-f_sim), NU * ALPHA),2);
    sim_XDS = mean(bsxfun(@times,bsxfun(@rdivide,f_sim,sum(f_sim+1)),(ALPHA+NU*SIGMA)),2);
end
