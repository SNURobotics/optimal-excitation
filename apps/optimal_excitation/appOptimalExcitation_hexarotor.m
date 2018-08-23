%% Trajectory Optimization for Dynamic Parameter Identification
close all
clear
clc

rmpath(genpath('apps/optimal_excitation/functions/'));
addpath(genpath('apps/optimal_excitation/functions_hexarotor/'));

%% Initialization
disp('initializing..')

robot     = makeHexarotor();           % robot model

trajectory.order           = 4;       % B Spline cubic base function
trajectory.horizon         = 0.1;      % trajectory horizon
trajectory.num_sample      = 10;     % number of samples of the trajectory
trajectory.base_frequency  = pi*0.3;  % Fourier trajectory base frequency

sigma     = eye(robot.dof);    % torque covariance
sigma_inv = pinv(sigma);

m         = 2;                       % num of trajectory coefs per joints
% p_initial = randc(m,6);
p_initial(:,1:3) = 0.01 *rand(m,3);  % initial p
p_initial(:,4:6) = 0.01*0.2*rand(m,3);
% % find feasible initial p
% while max(max(getConstraint(p_initial,trajectory,robot))) > 0
%     p_initial = 0*rand(m,robot.dof);
%     max(max(getConstraint(p_initial,trajectory,robot)))
% end
B_metric_inv_Phi_Bt_half = (robot.B_metric_inv_Phi_Bt)^(0.5);
B_metric_inv_Phi_Bt_half = (B_metric_inv_Phi_Bt_half+B_metric_inv_Phi_Bt_half')/2;
B_metric_inv_Phi_Bt_half_inv = pinv(B_metric_inv_Phi_Bt_half);
B_metric_inv_Phi_Bt_half_inv = (B_metric_inv_Phi_Bt_half_inv+B_metric_inv_Phi_Bt_half_inv')/2;
%% Trajectory Optimization

options = optimoptions(@fmincon,'Algorithm', 'sqp', 'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', false, 'TolCon',1e-7,'TolX',1e-5,'MaxFunEvals', ...
    1000000,'MaxIter',100000,'Display','iter','Hessian','bfgs'); %'fin-diff-grads'
% options = optimoptions(@fmincon,'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, 'Algorithm','sqp','TolCon',1e-15,'TolX',1e-15,'MaxFunEvals',1000000,'MaxIter',10000,'Display','iter');
% options = optimoptions(@fmincon, 'Algorithm','sqp','TolCon',1e-15,'TolX',1e-15,'MaxFunEvals',1000000,'MaxIter',10000,'Display','iter');
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp');

% w/ metric
disp('optimizing w/ metric..')

global C_previous
C_previous = zeros(10,10);
fval_previous = 1e10;
iter =1;
while(1)

    % [p_optimal,fval,exitflag,output,lam_costate] = fmincon(@(p)getCondNumber(p,robot,trajectory,sigma_inv), p_initial, [], [], [], [], [], [], @(p)getConstraint(p,trajectory,robot), options);
    [p_optimal,fval,exitflag, output,lam_costate] = fmincon(@(p)getCondNumber(p,robot,trajectory,sigma_inv), p_initial, [], [], [], [], [], [], [], options);
   
    p_initial(:,1:3) = 0.01 *rand(m,3);  % initial p
    p_initial(:,4:6) = 0.01*0.2*rand(m,3);

    if fval < fval_previous
        C_temp = getObjectiveMatrixC(p_optimal, robot, trajectory, sigma_inv);
        C_temp = B_metric_inv_Phi_Bt_half*C_temp*B_metric_inv_Phi_Bt_half_inv;
        C_temp = (C_temp+C_temp')/2;
        C_previous = (C_previous + C_temp);
        C_previous = C_previous * 0.2;
        fval_previous = fval;
        iter = iter+1
    end
    fval_previous
end
% w/o metric
robot.pd_metric_Phi = eye(10);
robot.B_metric_inv_Phi_Bt = eye(10);

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

for i = 1:num_calibration
    % Random trajecotry
    Phi_B_random_traj(:,i) = solveLeastSqaurePhiB(robot, p_initial, trajectory, sigma);
    
    % Optimization w/o metric
    Phi_B_wo_metric(:,i) = solveLeastSqaurePhiB(robot, p_baseline, trajectory, sigma);
    
    % Optimization w/  metric
    Phi_B_w_metric(:,i) = solveLeastSqaurePhiB(robot, p_optimal, trajectory, sigma);
end

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