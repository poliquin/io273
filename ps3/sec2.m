clear all
%%
rng(8675309);
dir = mfilename('fullpath');
cd(dir(1:end-4));

%% Simulate dataset
N=100;  K=2; SIGMA=[1,1;1,2]; ALPHA=1; BETA=1;
[Y, X, P,As] = sim_dataset(N, K, SIGMA, ALPHA, BETA);
disp('Mean and std. dev. for simulated data')
disp(table(num2str([mean(Y);std(Y)], '%3.3f &'), ...
    num2str([mean(X(:,1));std(X(:,1))], '%3.3f &'), ...
    num2str([mean(X(:,2));std(X(:,2))], '%3.3f &'), ...
            num2str([mean(P(:,1));std(P(:,1))], '%3.3f &'), ...
            num2str([mean(P(:,2));std(P(:,2))], '%3.3f \\\\'), ...
            'VariableNames', {'Y' 'X1' 'X2' 'P1' 'P2'}, ...
            'RowNames', {'Mean &' 'Std. Dev &'}))
disp(' ')

% Reshape X, P
X = reshape(X',[],1);
P = reshape(P',[],1);
%% Run gibbs using three different starting points and 10000 runs
runs = 5000; inits = 3;
init = zeros(inits,2); initsig = zeros(2,2,inits);
alpha = zeros(runs+1,inits); beta = zeros(runs+1,inits);
sigma = zeros(2,2,runs+1,inits);

% MCMC iterations
for ind = 1:inits
    init(ind,:) = mvnrnd([1,1],100*eye(2)); initsig(:,:,ind)= wishrnd(10*eye(2),100);
    printsig = initsig(:,:,ind)/initsig(1,1,ind);
    fprintf('Initial values: [beta, alpha, sig12, sig21, sig22] = \n %3.3f & %3.3f & %3.3f & %3.3f & %3.3f \\\\ \n',init(ind,1),init(ind,2),printsig(1,2),printsig(2,1),printsig(2,2))
    [A,B,S,YSTR] = gibbs(Y,[X,P],runs,init(ind,:),initsig(:,:,ind));
    scale = squeeze(S(1,1,:));
    beta(:,ind) = B./sqrt(scale);
    alpha(:,ind) = A./sqrt(scale);
    sigma(:,:,:,ind) = S./repmat(S(1,1,:),2);
end
%     beta = B./sqrt(scale);
%     alpha = A./sqrt(scale);
%     sigma = S./repmat(S(1,1,:),2);

% get confidence intervals
cut = (runs-1000)*0.025;
a_min=zeros(inits,1); a_max=zeros(inits,1); b_min=zeros(inits,1); b_max=zeros(inits,1);
s12_min=zeros(inits,1);s21_min=zeros(inits,1);s22_min=zeros(inits,1);
s12_max=zeros(inits,1);s21_max=zeros(inits,1);s22_max=zeros(inits,1);
for ind = 1:inits
    alpha_s = sort(alpha(3000:end,ind));
    a_min(ind) = min(alpha_s(cut:end-cut));
    a_max(ind) = max(alpha_s(cut:end-cut));
    beta_s = sort(beta(3000:end,ind));
    b_min(ind) = min(beta_s(cut:end-cut));
    b_max(ind) = max(beta_s(cut:end-cut));
    s12 = sort(sigma(1,2,3000:end,ind));
    s21 = sort(sigma(2,1,3000:end,ind));
    s22 = sort(sigma(2,2,3000:end,ind));
    s12_min(ind) = min(squeeze(s12(cut:end-cut)));
    s21_min(ind) = min(squeeze(s21(cut:end-cut)));
    s22_min(ind) = min(squeeze(s22(cut:end-cut)));
    s12_max(ind) = max(squeeze(s12(cut:end-cut)));
    s21_max(ind) = max(squeeze(s21(cut:end-cut)));
    s22_max(ind) = max(squeeze(s22(cut:end-cut)));
end

% display confidence intervals
for ind = 1:inits
    disp(fprintf('Mean and 95%% confidence interval for simulated coefficents, %i of %i',ind,inits))
    disp(table(...
        num2str([mean(alpha(3000:end,ind));a_min(ind);a_max(ind)], '%3.3f &'), ...
        num2str([mean(beta(3000:end,ind));b_min(ind);b_max(ind)], '%3.3f &'), ...
        num2str([mean(sigma(1,2,3000:end,ind));s12_min(ind);s12_max(ind)], '%3.3f &'), ...
        num2str([mean(sigma(2,2,3000:end,ind));s22_min(ind);s22_max(ind)], '%3.3f \\\\'), ...
        'VariableNames', {'Alpha' 'Beta' 's12' 's22'}, ...
        'RowNames', {'Mean &' 'lower 95%% &' 'upper 95%%'}))
    disp(' ')
end
%% Graphs


%% Part 2
clear all
[y2, x2, p2, z2] = sim_dataset2(100,2,1,1);
disp('Mean and std. dev. for simulated data (Part 2)')
disp(table(num2str([mean(y2);std(y2)], '%3.3f &'), ...
    num2str([mean(x2(:,1));std(x2(:,1))], '%3.3f &'), ...
    num2str([mean(x2(:,2));std(x2(:,2))], '%3.3f &'), ...
    num2str([mean(p2(:,1));std(p2(:,1))], '%3.3f &'), ...
    num2str([mean(p2(:,2));std(p2(:,2))], '%3.3f &'), ...
    num2str([mean(z2(:,1));std(z2(:,1))], '%3.3f &'), ...
    num2str([mean(z2(:,2));std(z2(:,2))], '%3.3f \\\\'), ...
    'VariableNames', {'Y' 'X1' 'X2' 'P1' 'P2' 'Z1' 'Z2'}, ...
    'RowNames', {'Mean &' 'Std. Dev &'}))
disp(' ')
y2_in = reshape(y2',[],1);
x2_in = reshape(x2',[],1);
p2_in = reshape(p2',[],1);
%% Part 2-1
inits = 3; runs =10000;
beta1=zeros(runs+1,inits); alpha1=zeros(runs+1,inits);
sigma1=zeros(2,2,runs+1,inits);
for ind = 1:inits
    init(ind,:) = mvnrnd([1,1],100*eye(2)); initsig(:,:,ind)= wishrnd(10*eye(2),100);
    printsig = initsig(:,:,ind)/initsig(1,1,ind);
    fprintf('Initial values: [beta, alpha, sig12, sig21, sig22] = \n %3.3f & %3.3f & %3.3f & %3.3f & %3.3f \\\\ \n',init(ind,1),init(ind,2),printsig(1,2),printsig(2,1),printsig(2,2))
    [a1,b1,s1,ystr] = gibbs(y2_in,[x2_in,p2_in],runs,init(ind,:),initsig(:,:,ind));
    scale = squeeze(s1(1,1,:));
    beta1(:,ind) = b1./sqrt(scale);
    alpha1(:,ind) = a1./sqrt(scale);
    sigma1(:,:,:,ind) = s1./repmat(s1(1,1,:),2);
end

% get confidence intervals
cut = (runs-1000)*0.025;
for ind = 1:inits
    alpha_s = sort(alpha1(1000:end,ind));
    a_min(ind) = min(alpha_s(cut:end-cut));
    a_max(ind) = max(alpha_s(cut:end-cut));
    beta_s = sort(beta1(1000:end,ind));
    b_min(ind) = min(beta_s(cut:end-cut));
    b_max(ind) = max(beta_s(cut:end-cut));
    s12 = sort(sigma1(1,2,:,ind));
    s21 = sort(sigma1(2,1,:,ind));
    s22 = sort(sigma1(2,2,:,ind));
    s12_min(ind) = min(squeeze(s12(cut:end-cut)));
    s21_min(ind) = min(squeeze(s21(cut:end-cut)));
    s22_min(ind) = min(squeeze(s22(cut:end-cut)));
    s12_max(ind) = max(squeeze(s12(cut:end-cut)));
    s21_max(ind) = max(squeeze(s21(cut:end-cut)));
    s22_max(ind) = max(squeeze(s22(cut:end-cut)));
end

% display confidence intervals
for ind = 1:inits
    disp(fprintf('Mean and 95%% confidence interval for simulated coefficents, %i of %i',ind,inits))
    disp(table(...
        num2str([mean(alpha(1000:end,ind));a_min(ind);a_max(ind)], '%3.3f &'), ...
        num2str([mean(beta(1000:end,ind));b_min(ind);b_max(ind)], '%3.3f &'), ...
        num2str([mean(sigma(1,2,1000:end,ind));s12_min(ind);s12_max(ind)], '%3.3f &'), ...
        num2str([mean(sigma(2,2,1000:end,ind));s22_min(ind);s22_max(ind)], '%3.3f \\\\'), ...
        'VariableNames', {'Alpha' 'Beta' 's12' 's22'}, ...
        'RowNames', {'Mean &' 'lower 95%% &' 'upper 95%%'}))
    disp(' ')
end
%% Part 2-2
% runs = 10000; inits=3;
% beta2=zeros(runs+1,inits); gamm2=zeros(runs+1,inits);
% sigma2=zeros(4,4,runs+1,inits);
% % Get dataset into desired form
% inputX = [reshape(x2',[],1),reshape(z2',[],1)];
% inputY = [y2, p2];
% for ind = 1:inits
%     initsig = zeros(4,4,inits);
%     init(ind,:) = mvnrnd([1,1],100*eye(2)); initsig(:,:,ind)= wishrnd(10*eye(4),100);
%     printsig = initsig(:,:,ind)/initsig(1,1,ind);
%     fprintf('Initial values: [beta, gamma] =\n %3.3f & %3.3f \\\\ \n',init(ind,1),init(ind,2))
%     [g2,b2,s2,ystr]=gibbs(inputY,inputX,runs,init(ind,:),initsig(:,:,ind));
%     scale2 = squeeze(s2(1,1,:));
%     beta2(:,ind) = b2./sqrt(scale2);
%     gamm2(:,ind) = g2./sqrt(scale2);
%     sigma2(:,ind) = s2./repmat(s2(1,1,:),4);
% end
% mean(b2(1000:end))
% mean(g2(1000:end))
% mean(s2(:,:,1000:end),3)


% GAMMA = 
