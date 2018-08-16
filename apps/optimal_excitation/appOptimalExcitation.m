%% Trajectory Optimization for Parameter Identification
close all
clear
clc

%% Initialization
robot     = makeWam7();         % robot model
m         = 10;                 % num of spline coefs per joints
initial_p = 10*rand(m,robot.dof);  % initial p

trajectory.order      = 4;      % cubic spline
trajectory.horizon    = 30;     % trajectory horizon
trajectory.num_sample = 100;    % number of samples of the trajectory

sigma = eye(robot.dof);         % torque covariance
sigma_inv = pinv(sigma);

%% Trajectory Optimization
% w/ metric
options = optimoptions(@fmincon,'Algorithm', 'interior-point', 'SpecifyObjectiveGradient', false, 'SpecifyConstraintGradient', false, 'TolCon',1e-7,'TolX',1e-5,'MaxFunEvals', ...
    1000000,'MaxIter',1000,'Display','iter'); %'fin-diff-grads'
% options = optimoptions(@fmincon,'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, 'Algorithm','sqp','TolCon',1e-15,'TolX',1e-15,'MaxFunEvals',1000000,'MaxIter',10000,'Display','iter');
% options = optimoptions(@fmincon, 'Algorithm','sqp','TolCon',1e-15,'TolX',1e-15,'MaxFunEvals',1000000,'MaxIter',10000,'Display','iter');
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
tic
[x,fval,exitflag,output,lam_costate] = fmincon(@(p)getCondNumber(p,robot,trajectory,sigma_inv), initial_p, [], [], [], [], [], [], [], options);
toc

% w/o metric
tic
[x_baseline,fvalx_baseline,exitflagx_baseline,outputx_baseline,lam_costatex_baseline] = fmincon(@(p)getCondNumberNoMetric(p,robot,trajectory,sigma_inv), initial_p, [], [], [], [], [], [], [], options);
toc

%% Optimization Result
initial_f = getCondNumber(initial_p, robot, trajectory,sigma_inv)
final_f   = getCondNumber(x, robot, trajectory, sigma_inv)
optimal_p = x

%% Least-Square Parameter Identification
% True Phi_B
Phi_B = robot.B * robot.Phi;

% Random trajecotry
Phi_B_random_traj = solveLeastSqaurePhiB(robot, initial_p, trajectory, sigma);

% Optimization w/o metric 
Phi_B_wo_metric = solveLeastSqaurePhiB(robot, x_baseline, trajectory, sigma);

% Optimization w/  metric
Phi_B_w_metric = solveLeastSqaurePhiB(robot, optimal_p, trajectory, sigma);

X = [Phi_B Phi_B_random_traj  Phi_B_wo_metric Phi_B_w_metric];