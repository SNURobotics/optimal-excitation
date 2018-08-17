%% Trajectory Optimization for Parameter Identification
close all
clear
clc

%% Initialization
robot     = makeWam7();         % robot model
m         = 10;                 % num of spline coefs per joints
p_initial = 20*rand(m,robot.dof) - 10;  % initial p

for i = 1:robot.dof
    p_initial(:,i) = 10*(rand(m,1)-0.5)/robot.Phi(10*(i-1)+1);
end

trajectory.order      = 4;      % cubic spline
trajectory.horizon    = 30;     % trajectory horizon
trajectory.num_sample = 100;    % number of samples of the trajectory

sigma = eye(robot.dof);         % torque covariance
sigma_inv = pinv(sigma);

%% Trajectory Optimization
% w/ metric
options = optimoptions(@fmincon,'Algorithm', 'interior-point', 'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', false, 'TolCon',1e-7,'TolX',1e-5,'MaxFunEvals', ...
    1000000,'MaxIter',1000,'Display','iter','Hessian','bfgs'); %'fin-diff-grads'
% options = optimoptions(@fmincon,'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, 'Algorithm','sqp','TolCon',1e-15,'TolX',1e-15,'MaxFunEvals',1000000,'MaxIter',10000,'Display','iter');
% options = optimoptions(@fmincon, 'Algorithm','sqp','TolCon',1e-15,'TolX',1e-15,'MaxFunEvals',1000000,'MaxIter',10000,'Display','iter');
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
tic
[p_optimal,fval,exitflag,output,lam_costate] = fmincon(@(p)getCondNumber(p,robot,trajectory,sigma_inv), p_initial, [], [], [], [], [], [], [], options);
toc

% w/o metric
% tic
% [p_baseline,fval_baseline,exitflag_baseline,output_baseline,lam_costate_baseline] = fmincon(@(p)getCondNumberNoMetric(p,robot,trajectory,sigma_inv), p_initial, [], [], [], [], [], [], [], options);
% toc
% 
% %% Least-Square Parameter Identification
% % True Phi_B
% Phi_B = robot.B * robot.Phi;
% 
% % Random trajecotry
% Phi_B_random_traj = solveLeastSqaurePhiB(robot, p_initial, trajectory, sigma);
% 
% % Optimization w/o metric 
% Phi_B_wo_metric = solveLeastSqaurePhiB(robot, p_baseline, trajectory, sigma);
% 
% % Optimization w/  metric
% Phi_B_w_metric = solveLeastSqaurePhiB(robot, p_optimal, trajectory, sigma);
% 
% X = [Phi_B Phi_B_random_traj  Phi_B_wo_metric Phi_B_w_metric];