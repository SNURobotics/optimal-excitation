%% Trajectory Optimization for Dynamic Parameter Identification
close all
clear
clc

rmpath(genpath('apps/optimal_excitation/functions_hexarotor/'));
rmpath(genpath('apps/optimal_excitation/functions/'));
addpath(genpath('apps/optimal_excitation/functions_Atlas/'));

%% Initialization
disp('initializing..')

robot     = makeKukaR820();           % robot model

trajectory.order           = 4;       % B Spline cubic base function
trajectory.horizon         = 20;      % trajectory horizon
trajectory.num_sample      = 1000;     % number of samples of the trajectory
trajectory.base_frequency  = pi*0.1;  % Fourier trajectory base frequency

sigma     = eye(robot.dof) * 1e-2;    % torque covariance
sigma_inv = pinv(sigma);

m         = 20;                       % num of trajectory coefs per joints
p_initial = rand(m,robot.dof) - 0.5;  % initial p

% find feasible initial p
while max(max(getConstraint(p_initial,trajectory,robot))) > 0
    p_initial = (rand(m,robot.dof) - 0.5)/5;
end

%% Trajectory Optimization

options = optimoptions(@fmincon,'Algorithm', 'sqp', 'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', false, 'TolCon',1e-7,'TolX',1e-5,'MaxFunEvals', ...
    1000000,'MaxIter',1000,'Display','iter','Hessian','bfgs'); %'fin-diff-grads'
% options = optimoptions(@fmincon,'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, 'Algorithm','sqp','TolCon',1e-15,'TolX',1e-15,'MaxFunEvals',1000000,'MaxIter',10000,'Display','iter');
% options = optimoptions(@fmincon, 'Algorithm','sqp','TolCon',1e-15,'TolX',1e-15,'MaxFunEvals',1000000,'MaxIter',10000,'Display','iter');
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp');

% w/ metric
disp('optimizing w/ metric..')
[p_optimal,fval,exitflag,output,lam_costate] = fmincon(@(p)getCondNumber(p,robot,trajectory,sigma_inv), p_initial, [], [], [], [], [], [], @(p)getConstraint(p,trajectory,robot), options);
% [p_optimal,fval,exitflag,output,lam_costate] = fmincon(@(p)getCondNumber(p,robot,trajectory,sigma_inv), p_initial, [], [], [], [], [], [], [], options);

% w/o metric
robot.pd_metric_Phi = eye(robot.dof*10);
robot.B_metric_inv_Phi_Bt = eye(39);

disp('optimizing w/o metric..')
[p_baseline,fval_baseline,exitflag_baseline,output_baseline,lam_costate_baseline] = fmincon(@(p)getCondNumber(p,robot,trajectory,sigma_inv), p_initial, [], [], [], [], [], [], @(p)getConstraint(p,trajectory,robot), options);


%% Least-Square Parameter Identification
disp('solving least square..')

% True Phi_B
Phi_B = robot.B * robot.Phi;
num_phi_b = size(Phi_B,1);

num_calibration = 100;

Phi_B_random_traj = zeros(num_phi_b,num_calibration);
Phi_B_wo_metric   = zeros(num_phi_b,num_calibration);
Phi_B_w_metric    = zeros(num_phi_b,num_calibration);

parfor i = 1:num_calibration
    % Random trajecotry
    Phi_B_random_traj(:,i) = solveLeastSqaurePhiB(robot, p_initial, trajectory, sigma);

    % Optimization w/o metric 
    Phi_B_wo_metric(:,i) = solveLeastSqaurePhiB(robot, p_baseline, trajectory, sigma);

    % Optimization w/  metric
    Phi_B_w_metric(:,i) = solveLeastSqaurePhiB(robot, p_optimal, trajectory, sigma);
    i
end

weighted_rms_error_random = sqrt(trace(pinv(robot.B_metric_inv_Phi_Bt)*(Phi_B_random_traj-robot.B*robot.Phi)*(Phi_B_random_traj-robot.B*robot.Phi)')/num_calibration);
weighted_rms_error_baseline = sqrt(trace(pinv(robot.B_metric_inv_Phi_Bt)*(Phi_B_wo_metric-robot.B*robot.Phi)*(Phi_B_wo_metric-robot.B*robot.Phi)')/num_calibration);
weighted_rms_error_optimal = sqrt(trace(pinv(robot.B_metric_inv_Phi_Bt)*(Phi_B_w_metric-robot.B*robot.Phi)*(Phi_B_w_metric-robot.B*robot.Phi)')/num_calibration);


variance_random = var(Phi_B_random_traj, 0 ,2);
variance_baseline = var(Phi_B_wo_metric, 0 ,2);
variance_optimal = var(Phi_B_w_metric, 0, 2);
variance = [variance_random variance_baseline variance_optimal];
variance_proportion = [variance(:,1)./variance(:,1) variance(:,2)./variance(:,1) variance(:,3)./variance(:,1)];

% X = [Phi_B Phi_B_random_traj  Phi_B_wo_metric Phi_B_w_metric];
disp('..done')

%% Trajectory Plot
time_step = 0.1;
num_time = floor(trajectory.horizon / time_step);
sample_time = linspace(0,trajectory.horizon,num_time);

q_initial = makeFourier(p_initial, trajectory.base_frequency, sample_time);
q_optimal = makeFourier(p_optimal, trajectory.base_frequency, sample_time);

figure(); hold on;
for i = 1:robot.dof
    subplot(1,robot.dof,i);
    plot(q_initial(i,:),'b'); hold on;
    plot(q_optimal(i,:),'r');    
end

% robot visualization
appVisualizeKUKA(q_optimal)