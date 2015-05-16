% Problem Set 4
% =============
% Do Yoon Kim, Chris Poliquin, David Zhang
% May 15, 2015

rng(8675309);
dir = mfilename('fullpath');
cd(dir(1:end-4));
options = optimset('Display', 'final', 'TolFun', 10e-10, ...
                   'MaxIter', 1000, 'MaxFunEvals', 1000);


%% Section 2, Question 3
% ----------------------------------------------------------------------------
% Define true parameters
beta = .99;
RC = 10;
theta1  = .05;
theta30 = .3;
theta31 = .5;
theta32 = .2;

% Calculate expected continuation value
EV = ev(beta, RC, theta1, theta30, theta31, theta32);


%% Section 2, Question 4
% ----------------------------------------------------------------------------
% Plot continuation value

EVplot = EV(1:11,:);
EVplotx = linspace(0,10,11);
f = figure('PaperPosition', [.1, .2, 6.2, 3.5], 'PaperSize', [6.4, 4]);
plot(EVplotx, EVplot(:,1)', '--', EVplotx, EVplot(:, 2), ':')
title('Plot of Continuation Values EV')
xlabel('x')
ylabel('EV')
legend('i = 0', 'i = 1')
saveas(f, 'figs/EV.pdf');


%% Section 2, Question 5
% ----------------------------------------------------------------------------
% Simulate data

[xt, it] = sim_data(EV, beta, RC, theta1, theta30, theta31);
% Generate summary statistics and plots
xtit = zeros(100000, 2);
count = 1;
for i=1:size(xt, 1)
    for j=1:size(xt, 2)
        xtit(count, 1) = xt(i, j);
        xtit(count, 2) = it(i, j);
        count = count + 1;
    end
end

f = figure('PaperPosition', [.1, .2, 6.2, 3.5], 'PaperSize', [6.4, 4]);
histogram(xtit(xtit(:, 2)==1, 1))
title('Histogram of x when i=0')
xlabel('x')
ylabel('freq')
saveas(f, 'figs/histx.pdf');

for i=1:max(xtit(xtit(:, 2)==1,1))
    frac_i(i) = sum(xtit(xtit(:,2)==1 & xtit(:,1)==i,1))/sum(xtit(xtit(:,1)==i,1));
end

f = figure('PaperPosition', [.1, .2, 6.2, 3.5], 'PaperSize', [6.4, 4]);
bar(frac_i)
title('Probability of Investment')
xlabel('x')
ylabel('Fraction with i=1')
saveas(f, 'figs/prob.pdf');

sprintf('2.5\nMin = %f, \n25th Pctile = %f, \nMean = %f, \nMedian = %f, \n75th Pctile = %f, \nMax = %f',...
    min(xtit(xtit(:,2)==1,1)), prctile(xtit(xtit(:,2)==1,1),25),mean(xtit(xtit(:,2)==1,1)),median(xtit(xtit(:,2)==1,1)), prctile(xtit(xtit(:,2)==1,1),75), max(xtit(xtit(:,2)==1,1)))


%% Section 3, Question 1
% ----------------------------------------------------------------------------
% columnize choices and transitions
renew = reshape(it(1:999, :), [], 1);
trans = reshape(diff(xt), [], 1);

% estimate transition probabilities using cases without renewal
theta3 = tabulate(trans(renew == 0));
theta3 = theta3(:, 3) / 100;
theta3se = sqrt(theta3.*(1-theta3) ./ sum(renew == 0));

%% Use nested fixed point to find theta1
obj = @(theta1hat) -1 * rust(theta1hat, theta3(2), theta3(3), RC, beta, xt, it);
[theta1hat, ~, ~, ~, ~, hess] = fminunc(obj, rand, options);
se = sqrt(diag(inv(hess)));

sprintf('Theta 3')
sprintf('%0.04f (%0.04f)\n', [theta3, theta3se]')
sprintf('Theta 1 Rust:  %0.04f (%0.04f)', theta1hat, se)


%% Section 3, Question 2
% ----------------------------------------------------------------------------
% estimate theta_1 using conditional choice probabilities
theta1 = hotz(xt, it, RC, theta3, beta, options);

% bootstrap the standard errors
boots = bootstrp(500, @(x, i) hotz(x, i, RC, theta3, beta, options), xt, it);
sprintf('Theta 1 CCP:  %0.04f (%0.04f)', theta1, std(boots))

f = figure('PaperPosition', [.1, .2, 6.2, 3.5], 'PaperSize', [6.4, 4]);
histogram(boots, 20)
title('Estimates of $\theta_1$, 500 Replications', 'interpreter', 'latex')
xlabel('$\theta_1$', 'interpreter', 'latex')
saveas(f, 'figs/theta_1_hist.pdf');


%% Section 3, Question 3
% ----------------------------------------------------------------------------
EVNew = ev(beta, 20, .02, 1-theta31hat-theta32hat, theta31hat, theta32hat);
EVNew(1,1)
EVOld = ev(beta, 10, theta1hat, 1-theta31hat-theta32hat, theta31hat, theta32hat);
EVOld(1,1)


%% Section 3, Question 4
% ----------------------------------------------------------------------------
% Simulate replacement decisions for new engines

ERNew=zeros(30,1);
for i=1:30
    EVNew = ev(beta, i, .02, 1-theta31hat-theta32hat, theta31hat, theta32hat);
    pi=fppi(EVNew,beta, i, .02, 1-theta31hat-theta32hat, theta31hat, theta32hat);
    ERNew(i) = 100*sum(pi(:,2));
end

% Simulate replacement decisions for old engines. 
EROld=zeros(30,1);
for i=1:30
    EVOld = ev(beta, i, theta1hat, 1-theta31hat-theta32hat, theta31hat, theta32hat);
    pi=fppi(EVOld,beta, i, theta1hat, 1-theta31hat-theta32hat, theta31hat, theta32hat);
    EROld(i) = 100*sum(pi(:,2));
end

% Define xaxis
xaxis = linspace(1,30,30);
f = figure('PaperPosition', [.1, .2, 6.2, 3.5], 'PaperSize', [6.4, 4]);
plot(ERNew,xaxis','--',EROld,xaxis,':')
axis([1 30 0 32])
title('Plot of demand functions for new and old engines')
xlabel('RC')
ylabel('Expected number of replacements per period')
legend('New engine', 'Old engine')
saveas(f, 'figs/er.pdf');


%% Section 3, Question 5
% ----------------------------------------------------------------------------
% Simulate replacement decisions for new engines
EVOld = ev(beta, 10, theta1hat, 1-theta31hat-theta32hat, theta31hat, theta32hat);

EVDiff=zeros(30,1);
for i=1:30
    EVNew = ev(beta, i, .02, 1-theta31hat-theta32hat, theta31hat, theta32hat);
    EVDiff(i)=(EVNew(1,1)-EVOld(1,1))*100;
end

f = figure('PaperPosition', [.1, .2, 6.2, 3.5], 'PaperSize', [6.4, 4]);
plot(xaxis,EVDiff)
title('Total value difference between old and new engines for 100 buses')
xlabel('RC(new)')
ylabel('Total value difference')
saveas(f, 'figs/evdiff.pdf');
