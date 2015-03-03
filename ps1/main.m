% 14.273  Problem Set 1
% ----------------------------------------------------------------------------
% Do Yoon Kim & Christopher Poliquin

%% Section 1, Generate data
% ----------------------------------------------------------------------------
disp('Creating three datasets')
run('create_data.m')

%% Section 2, Demand-side estimation
% ----------------------------------------------------------------------------
clear all
disp('Section 2, Question 1')
run('sec2q1.m')  % Question 1
disp('Section 2, Question 2')
run('sec2q2.m')  % Question 2

%% Section 3, Adding supply-side moments
% ----------------------------------------------------------------------------
clear all
disp('Section 3, Question 1')
run('sec3q1.m')   % Question 1
disp('Section 3, Question 2')
run('sec3q2.m')   % Question 2, estimate demand and supply side
disp('Section 3, Question 2(b)')
run('sec3q2b.m')  % Question 2(b)
disp('Section 3, Question 2(c)')
sec3q2c();        % Question 2(c)

% Section 4, Merger simulation
% ----------------------------------------------------------------------------
clear all
disp('Section 4, Merger Simulation')
run('sec4.m')

