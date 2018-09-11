%% Trajectory Optimization for Dynamic Parameter Identification
close all
clear
clc

warning_handler = warning('query','last');
rmpath(genpath('apps/optimal_excitation/functions_KUKA/'));
rmpath(genpath('apps/optimal_excitation/functions_Atlas/'));
addpath(genpath('apps/optimal_excitation/functions_hexarotor/'));
warning('off',warning_handler.identifier);

%% Initialization
disp('initializing..')

robot     = makeHexarotor();           % robot model

trajectory.order           = 4;       % B Spline cubic base function
trajectory.horizon         = 5;      % trajectory horizon
trajectory.num_sample      = 10;     % number of samples of the trajectory
trajectory.base_frequency  = pi*0.3;  % Fourier trajectory base frequency

sigma     = eye(robot.dof);    % torque covariance
sigma_inv = pinv(sigma);

m         = 3;                       % num of trajectory coefs per joints
p_initial = rand(3,6);
% p_initial(2,1:3) = (rand(1,3)-0.5)*6;
% p_initial(2,4:6) = (rand(1,3)-0.5)*8*pi;
% p_initial(3,1:3) = (rand(1,3)-0.5)*2;
% p_initial(3,4:6) = (rand(1,3)-0.5)*4;

% B_metric_inv_Phi_Bt_half = (robot.B_metric_inv_Phi_Bt)^(0.5);
% B_metric_inv_Phi_Bt_half = (B_metric_inv_Phi_Bt_half+B_metric_inv_Phi_Bt_half')/2;
% B_metric_inv_Phi_Bt_half_inv = pinv(B_metric_inv_Phi_Bt_half);
% B_metric_inv_Phi_Bt_half_inv = (B_metric_inv_Phi_Bt_half_inv+B_metric_inv_Phi_Bt_half_inv')/2;


%% Trajectory Optimization

options = optimoptions(@fmincon,'Algorithm', 'sqp', 'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', false, 'TolCon',1e-7,'TolX',1e-5,'MaxFunEvals', ...
    1000000,'MaxIter',100000,'Display','off','Hessian','bfgs'); %'fin-diff-grads'
% options = optimoptions(@fmincon,'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, 'Algorithm','sqp','TolCon',1e-15,'TolX',1e-15,'MaxFunEvals',1000000,'MaxIter',10000,'Display','iter');
% options = optimoptions(@fmincon, 'Algorithm','sqp','TolCon',1e-15,'TolX',1e-15,'MaxFunEvals',1000000,'MaxIter',10000,'Display','iter');
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp');

% w/ metric2
disp('optimizing w/ metric..')

C_previous = zeros(10,10);
fval_previous = 1e20;
num_trajectory = 0;
discount = 1;

for iter = 1:1000000
    p_initial = rand(3,6);
%     p_initial(2,1:3) = (rand(1,3)-0.5)*6;
%     p_initial(2,4:6) = (rand(1,3)-0.5)*8*pi;
%     p_initial(3,1:3) = (rand(1,3)-0.5)*2;
%     p_initial(3,4:6) = (rand(1,3)-0.5)*4;   
    
    [p_optimal,fval,exitflag, output,lam_costate] = fmincon(@(p)getCondNumber(p,robot,trajectory,sigma_inv,C_previous), p_initial, [], [], [], [], [], [], @(p)getConstraint(p,trajectory,robot), options);

    if fval < fval_previous && exitflag ~= -2
        num_trajectory = num_trajectory+1;
        p_optimal_cummulative(:,:,num_trajectory) = p_optimal;
        
        C_temp = getObjectiveMatrixC(p_optimal, robot, trajectory, sigma_inv);
        C_temp = (C_temp+C_temp')/2;
        C_previous = (C_previous + C_temp)*discount;
        fval_previous = fval;
    end
    disp(['iter: ' num2str(iter) ' num_trajectory: ' num2str(num_trajectory) ' fval ' num2str(fval_previous)]);
end

p_optimal_cummulative(1,4:6,:) = rem(p_optimal_cummulative(1,4:6,:),2*pi);
p_optimal_cummulative(1,1:3,:) = 0;

C_previous = zeros(10,10);
for i = 1:num_trajectory-1
    C_temp = getObjectiveMatrixC(p_optimal_cummulative(:,:,i), robot, trajectory, sigma_inv);
    C_temp = (C_temp+C_temp')/2;
    C_previous = (C_previous + C_temp)*discount;
end
cond_number = getCondNumber(p_optimal_cummulative(:,:,num_trajectory),robot,trajectory,sigma_inv,C_previous)

% w/o metric
% robot.pd_metric_Phi = eye(10);
% robot.B_metric_inv_Phi_Bt = eye(10);
% 
% disp('optimizing w/o metric..')
% [p_baseline,fval_baseline,exitflag_baseline,output_baseline,lam_costate_baseline] = fmincon(@(p)getCondNumber(p,robot,trajectory,sigma_inv), p_initial, [], [], [], [], [], [], @(p)getConstraint(p,trajectory,robot), options);


%% Least-Square Parameter Identification
% disp('solving least square..')
% 
% % True Phi_B
% Phi_B = robot.B * robot.Phi;
% num_phi_b = size(Phi_B,1);
% 
% num_calibration = 100;
% 
% Phi_B_random_traj = zeros(num_phi_b,num_calibration);
% Phi_B_wo_metric   = zeros(num_phi_b,num_calibration);
% Phi_B_w_metric    = zeros(num_phi_b,num_calibration);
% 
% for i = 1:num_calibration
%     % Random trajecotry
%     Phi_B_random_traj(:,i) = solveLeastSqaurePhiB(robot, p_initial, trajectory, sigma);
%     
%     % Optimization w/o metric
%     Phi_B_wo_metric(:,i) = solveLeastSqaurePhiB(robot, p_baseline, trajectory, sigma);
%     
%     % Optimization w/  metric
%     Phi_B_w_metric(:,i) = solveLeastSqaurePhiB(robot, p_optimal, trajectory, sigma);
% end
% 
% variance_random = var(Phi_B_random_traj, 0 ,2);
% variance_baseline = var(Phi_B_wo_metric, 0 ,2);
% variance_optimal = var(Phi_B_w_metric, 0, 2);
% variance = [variance_random variance_baseline variance_optimal];
% variance_proportion = [variance(:,1)./variance(:,1) variance(:,2)./variance(:,1) variance(:,3)./variance(:,1)];
% 
% % X = [Phi_B Phi_B_random_traj  Phi_B_wo_metric Phi_B_w_metric];
% disp('..done')

%% Trajectory Plot
time_step = 0.1;
num_time = floor(trajectory.horizon / time_step);
sample_time = linspace(0,trajectory.horizon,num_time);

optimal_trajectory = zeros(robot.dof, (num_trajectory-1)*num_time);
for i = 1:num_trajectory-1
    optimal_trajectory(:, (i-1)*num_time+1:i*num_time) = make5thSplineP2P(p_optimal_cummulative(1,:,i)',...
                                                                  p_optimal_cummulative(2,:,i)',...
                                                                  p_optimal_cummulative(3,:,i)',...
                                                                  p_optimal_cummulative(1,:,i+1)',...
                                                                  p_optimal_cummulative(2,:,i+1)',...
                                                                  p_optimal_cummulative(3,:,i+1)',...
                                                                  trajectory.horizon, sample_time);
end

figure(); hold on;
total_time = linspace(0,trajectory.horizon*(num_trajectory-1),num_time*(num_trajectory-1));
for i = 1:robot.dof
    subplot(1,robot.dof,i);
    plot(total_time, optimal_trajectory(i,:)); hold on;

    plot(0, p_optimal_cummulative(1,i,1), '*b');
    for j = 1:(num_trajectory-1)
        plot(trajectory.horizon*j, p_optimal_cummulative(1,i,j+1), '*b');
    end
    hold off;
end

%%
appVisualizeHexarotor(optimal_trajectory);