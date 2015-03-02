% BLP estimation using demand and cost moments

%% Load the (3, 100) dataset
% ----------------------------------------------------------------------------
load('data/100_3.mat')
prodcount = 3;
mktcount = 100;

% if dataset has zero shares, we need to replace the values
tabulate(shares == 0)
shares(shares == 0) = 10e-6;
prices(shares == 0) = 0;
prods(shares == 0, :) = 0;

%% Run the estimation procedure using cost and demand side moments
% ----------------------------------------------------------------------------
runs = 30;
% run the procedure multiple times so that we get multiple starting values
[estimates, variances, fval, true_etas, siml_etas] = blpdemand(prices, prods, shares, cost, ...
            prodcount, mktcount, true, true, runs);

% show all of the estimates
disp('Full results')
disp(estimates)

% find best set of coefficients and standard errors (smallest fval)
[fval, minidx] = min(estimates(:, 1));
theta = estimates(minidx, 2:end);
vcov  = reshape(variances(minidx, :), 5, 5);

% calculate bias of the estimates
bias = theta' - [-1; 5; 1; 1; 1];

% print the best estimates
disp(' ')
disp('(3, 100) Coefficients and Standard Errors (with Cost Shifter)')
disp(table(num2str(theta', '%3.3f &'), num2str(sqrt(diag(vcov)), '(%3.3f) &'), ...
            num2str(bias, '%4.3f \\\\'), ...
            'VariableNames', {'Coef' 'SE' 'Bias'}, ...
            'RowNames', {'Price  &', 'X1  &', 'X2  &', 'X3  &', 'Sigma  &'}))
fprintf('Objective function: %i', fval)
disp(' ')

% print the (average) elasticities 
true_eta = reshape(true_etas(minidx,:),prodcount,[]);
siml_eta = reshape(siml_etas(minidx,:),prodcount,[]);
final_true_eta = mean(true_eta,2);
final_siml_eta = mean(siml_eta,2);
disp(' ')
disp('(3, 100) Own Price Elasticities')
disp(table(num2str(final_true_eta, '%3.3f & '), ...
            num2str(final_siml_eta, '%3.3f \\\\'), ...
            'VariableNames',{'True', 'Simlulated'},...
            'RowNames',{'Firm 1 &','Firm 2 &','Firm 3 &'}))

% print elasticities across markets
disp(table(num2str([true_eta;siml_eta], '%3.3f & '), ...
            'RowNames',{'Firm 1 True &','Firm 2 True &','Firm 3 True &',...
            'Firm 1 Siml &','Firm 2 Siml &','Firm 3 Siml &'}))
        
% Average profits
true_pi = mean(-reshape(prices,3,[]).*reshape(shares,3,[])./true_eta*500,2);
siml_pi = mean(-reshape(prices,3,[]).*reshape(shares,3,[])./siml_eta*500,2);

% Average surplus
true_cs = mean(reshape(([prices, prods,-prices*0.5] * [-1;5;1;1;1]).*shares,3,[]),2)
siml_cs = mean(reshape(([prices, prods,prices*theta(1)*0.5] * theta').*shares,3,[]),2)

disp(table(num2str(true_pi,'%4.3f &'),num2str(siml_pi,'%3.3f &'),...
    num2str(true_cs,'%4.3f &'), num2str(siml_cs, '%3.3f \\\\ '),...
    'VariableNames',{'True_profit','Siml_profit','True_CS','Siml_CS'},...
    'RowNames',{'Firm 1 &','Firm 2 &','Firm 3 &'}))

%% Load the (3, 10) dataset
% ----------------------------------------------------------------------------
load('data/10_3.mat')
prodcount = 3;
mktcount = 10;

% if dataset has zero shares, we need to replace the values
tabulate(shares == 0)
shares(shares == 0) = 10e-6;
prices(shares == 0) = 0;
prods(shares == 0, :) = 0;

%% Run the estimation procedure using cost and demand side moments
% ----------------------------------------------------------------------------
runs = 30;
% run the procedure multiple times so that we get multiple starting values
[estimates, variances, fval, true_etas, siml_etas] = blpdemand(prices, prods, shares, cost, ...
            prodcount, mktcount, true, true, runs);

% show all of the estimates
disp('Full results')
disp(estimates)

% find best set of coefficients and standard errors (smallest fval)
[fval, minidx] = min(estimates(:, 1));
theta = estimates(minidx, 2:end);
vcov  = reshape(variances(minidx, :), 5, 5);

% calculate bias of the estimates
bias = theta' - [-1; 5; 1; 1; 1];

% print the best estimates
disp(' ')
disp('(3, 10) Coefficients and Standard Errors (with Cost Shifter)')
disp(table(num2str(theta', '%3.3f &'), num2str(sqrt(diag(vcov)), '(%3.3f) &'), ...
            num2str(bias, '%4.3f \\\\'), ...
            'VariableNames', {'Coef' 'SE' 'Bias'}, ...
            'RowNames', {'Price  &', 'X1  &', 'X2  &', 'X3  &', 'Sigma  &'}))
fprintf('Objective function: %i', fval)
disp(' ')

% print the (average) elasticities 
true_eta = reshape(true_etas(minidx,:),prodcount,[]);
siml_eta = reshape(siml_etas(minidx,:),prodcount,[]);
final_true_eta = mean(true_eta,2);
final_siml_eta = mean(siml_eta,2);
disp(' ')
disp('(3, 100) Own Price Elasticities')
disp(table(num2str(final_true_eta, '%3.3f & '), ...
            num2str(final_siml_eta, '%3.3f \\\\'), ...
            'VariableNames',{'True', 'Simlulated'},...
            'RowNames',{'Firm 1 &','Firm 2 &','Firm 3 &'}))

% print elasticities across markets
disp(table(num2str([true_eta;siml_eta], '%3.3f & '), ...
            'RowNames',{'Firm 1 True &','Firm 2 True &','Firm 3 True &',...
            'Firm 1 Siml &','Firm 2 Siml &','Firm 3 Siml &'}))
        
% Average profits
true_pi = mean(-reshape(prices,3,[]).*reshape(shares,3,[])./true_eta*500,2);
siml_pi = mean(-reshape(prices,3,[]).*reshape(shares,3,[])./siml_eta*500,2);

% Average surplus
true_cs = mean(reshape(([prices, prods,-prices*0.5] * [-1;5;1;1;1]).*shares,3,[]),2)
siml_cs = mean(reshape(([prices, prods,prices*theta(1)*0.5] * theta').*shares,3,[]),2)

disp(table(num2str(true_pi,'%4.3f &'),num2str(siml_pi,'%3.3f &'),...
    num2str(true_cs,'%4.3f &'), num2str(siml_cs, '%3.3f \\\\ '),...
    'VariableNames',{'True_profit','Siml_profit','True_CS','Siml_CS'},...
    'RowNames',{'Firm 1 &','Firm 2 &','Firm 3 &'}))
