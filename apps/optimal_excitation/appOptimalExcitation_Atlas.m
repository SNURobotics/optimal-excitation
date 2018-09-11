%% Trajectory Optimization for Dynamic Parameter Identification
close all
clear
clc

warning_handler = warning('query','last');
rmpath(genpath('apps/optimal_excitation/functions_hexarotor/'));
rmpath(genpath('apps/optimal_excitation/functions_KUKA/'));
addpath(genpath('apps/optimal_excitation/functions_Atlas/'));
warning('off',warning_handler.identifier);

%% Initialization
disp('initializing..')

robot     = makeAtlasV5();           % robot model

trajectory.order           = 4;       % B Spline cubic base function
trajectory.horizon         = 20;      % trajectory horizon
trajectory.num_sample      = 50;     % number of samples of the trajectory
trajectory.base_frequency  = pi*0.1;  % Fourier trajectory base frequency

sigma     = eye(6) * 1e-2;    % torque covariance
sigma_inv = pinv(sigma);

m         = 10;                       % num of trajectory coefs per joints
p_initial = rand(m,robot.dof);  % initial p

% find feasible initial p
tries = 0;
while max(max(getConstraint(p_initial,trajectory,robot))) > 0 && getCondNumber(p_initial,robot,trajectory,sigma_inv) < 0
    tries = tries + 1;
    disp([num2str(tries) 'th initial guess']);
    p_initial = (rand(m,robot.dof))/30;
end

%% Trajectory Optimization

options = optimoptions(@fmincon,'Algorithm', 'sqp', 'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', false, 'TolCon',1e-7,'TolX',1e-5,'MaxFunEvals', ...
    1000000,'MaxIter',1000,'Display','iter','Hessian','bfgs'); %'fin-diff-grads'
% options = optimoptions(@fmincon,'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, 'Algorithm','sqp','TolCon',1e-15,'TolX',1e-15,'MaxFunEvals',1000000,'MaxIter',10000,'Display','iter');
% options = optimoptions(@fmincon, 'Algorithm','sqp','TolCon',1e-15,'TolX',1e-15,'MaxFunEvals',1000000,'MaxIter',10000,'Display','iter');
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
tic
% w/ metric
disp('optimizing w/ metric..')
[p_optimal,fval,exitflag,output,lam_costate] = fmincon(@(p)getCondNumber(p,robot,trajectory,sigma_inv), p_initial, [], [], [], [], [], [], @(p)getConstraint(p,trajectory,robot), options);
% [p_optimal,fval,exitflag,output,lam_costate] = fmincon(@(p)getCondNumber(p,robot,trajectory,sigma_inv), p_initial, [], [], [], [], [], [], [], options);
toc
% w/o metric
robot.pd_metric_Phi = eye(robot.dof*10);
robot.B_metric_inv_Phi_Bt = eye(201);

disp('optimizing w/o metric..')
[p_baseline,fval_baseline,exitflag_baseline,output_baseline,lam_costate_baseline] = fmincon(@(p)getCondNumber(p,robot,trajectory,sigma_inv), p_initial, [], [], [], [], [], [], @(p)getConstraint(p,trajectory,robot), options);
