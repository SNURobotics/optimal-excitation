%% Trajectory Optimization for Dynamic Parameter Identification: Find Only the Sufficiently Excited Parameters
close all
clear
clc

warning_handler = warning('query','last');
rmpath(genpath('apps/optimal_excitation/functions_hexarotor/'));
rmpath(genpath('apps/optimal_excitation/functions_Atlas/'));
addpath(genpath('apps/optimal_excitation/functions_KUKA/'));
warning('off',warning_handler.identifier);

%% Options
plot_and_visualize = false;
use_cond_number    = false;

%% Initialization
disp('initializing..')

robot     = makeKukaR820();           % robot model

trajectory.order           = 4;       % B Spline cubic base function
trajectory.horizon         = 10;      % trajectory horizon
trajectory.num_sample      = 1000;     % number of samples of the trajectory
trajectory.base_frequency  = pi*0.3;  % Fourier trajectory base frequency

sigma     = eye(robot.dof) * 1e-2;    % torque covariance
sigma_inv = inv(sigma);

max_fval = inf;

m         = 10;                       % num of trajectory coefs per joints
p_initial = rand(m,robot.dof) - 0.5;  % initial p

% find feasible initial p
while max(max(getConstraint(p_initial,trajectory,robot))) > 0
    p_initial = (rand(m,robot.dof) - 0.5)/10;
end
p_initial_0 = p_initial;

% nominal phi
B = robot.B;
num_base = size(B,1);

robot.Phi_0 = zeros(10*robot.dof,1);  % nominal phi
for i = 1:robot.dof
    robot.Phi_0(1+(i-1)*10:i*10) = getNominalPhi(robot.Phi(1+(i-1)*10:i*10));
end

robot.pd_metric_Phi_0 = zeros(10*robot.dof, 10*robot.dof);
for i = 1:robot.dof
    robot.pd_metric_Phi_0(10*(i-1)+1:10*i, 10*(i-1)+1:10*i) = getPDMetricInertiaPhi(robot.Phi_0(10*(i-1)+1:10*i));
end

robot.B_metric_inv_Phi_Bt_0 = B * pinv(robot.pd_metric_Phi_0) * B';
robot.B_metric_inv_Phi_Bt_0 = (robot.B_metric_inv_Phi_Bt_0 + robot.B_metric_inv_Phi_Bt_0')/2;


%% Trajectory Optimization

options = optimoptions(@fmincon,'Algorithm', 'sqp', 'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', false, 'TolCon',1e-7,'TolX',1e-5,'MaxFunEvals', ...
    1000000,'MaxIter',1000,'Display','iter','Hessian','bfgs'); %'fin-diff-grads'

% w/ metric
disp('optimizing w/ metric..')

num_opt = 1     % optimize [num_opt]th eigen value
while true
    if use_cond_number
        [p_optimal,fval,exitflag,output,lam_costate] = fmincon(@(p)getCondNumber(p,robot,trajectory,sigma_inv,num_opt), p_initial, [], [], [], [], [], [], @(p)getConstraint(p,trajectory,robot), options);
    else
        [p_optimal,fval,exitflag,output,lam_costate] = fmincon(@(p)getEOptimality(p,robot,trajectory,sigma_inv,num_opt), p_initial, [], [], [], [], [], [], @(p)getConstraint(p,trajectory,robot), options);
    end
    
    if fval < max_fval
        break;
    else
        num_opt = num_opt + 1
        p_initial = p_optimal;
    end
end

% w/o metric
disp('optimizing w/o metric..')

robot.B_metric_inv_Phi_Bt_0 = eye(size(robot.B_metric_inv_Phi_Bt));
p_initial = p_initial_0;

num_opt = 1     % optimize [num_opt]th eigen value
while true
    if use_cond_number
        [p_baseline,fval,exitflag,output,lam_costate] = fmincon(@(p)getCondNumber(p,robot,trajectory,sigma_inv,num_opt), p_initial, [], [], [], [], [], [], @(p)getConstraint(p,trajectory,robot), options);
    else
        [p_baseline,fval,exitflag,output,lam_costate] = fmincon(@(p)getEOptimality(p,robot,trajectory,sigma_inv,num_opt), p_initial, [], [], [], [], [], [], @(p)getConstraint(p,trajectory,robot), options);
    end
    
    if fval < max_fval
        break;
    else
        num_opt = num_opt + 1
        p_initial = p_optimal;
    end
end

%% Least-Square Parameter Identification
disp('solving least square..')

load('p_test_kuka.mat');
% p_test = (rand(m,robot.dof) - 0.5);


for iter = 1 : 100
% w/ metric
robot.B_metric_inv_Phi_Bt_0 = B * pinv(robot.pd_metric_Phi_0) * B';
robot.B_metric_inv_Phi_Bt_0 = (robot.B_metric_inv_Phi_Bt_0 + robot.B_metric_inv_Phi_Bt_0')/2;

C = getObjectiveMatrixC(p_optimal, robot, trajectory, sigma_inv);
[Q, D] = eig(C);
D = inv(D);

k = size(D,1);
eigenvalues = zeros(k,1);
for i=1:k
    eigenvalues(i) = D(i,i);    
end

[~, indicies] = sort(eigenvalues, 'descend');
index_opt = indicies(num_opt,1);
max_variance = eigenvalues(index_opt);

Phi_B_w_metric = solveLeastSqaurePhiB(robot,p_optimal,trajectory,sigma,max_variance);

% w/o metric
robot.B_metric_inv_Phi_Bt_0 = eye(size(robot.B_metric_inv_Phi_Bt));

C = getObjectiveMatrixC(p_baseline, robot, trajectory, sigma_inv);
[Q, D] = eig(C);
D = inv(D);

k = size(D,1);
eigenvalues = zeros(k,1);
for i=1:k
    eigenvalues(i) = D(i,i);    
end

[~, indicies] = sort(eigenvalues, 'descend');
index_opt = indicies(num_opt,1);
max_variance = eigenvalues(index_opt);

Phi_B_wo_metric = solveLeastSqaurePhiB(robot,p_baseline,trajectory,sigma,max_variance);


Phi_B_initial = solveLeastSqaurePhiB(robot, p_initial_0, trajectory, sigma,inf);

err_w_metric(iter) = norm((robot.B_metric_inv_Phi_Bt)^(-0.5)*(Phi_B_w_metric-B*robot.Phi))./sqrt(size(B,1));
err_wo_metric(iter) = norm((robot.B_metric_inv_Phi_Bt)^(-0.5)*(Phi_B_wo_metric-B*robot.Phi))./sqrt(size(B,1));
err_initial(iter) = norm((robot.B_metric_inv_Phi_Bt)^(-0.5)*(Phi_B_initial-B*robot.Phi))./sqrt(size(B,1));

[A_B, b] = getFullRegressorTorque(robot, p_test, trajectory);
torque_err_w_metric(iter) = rms(A_B*Phi_B_w_metric - b);
torque_err_wo_metric(iter) = rms(A_B*Phi_B_wo_metric - b);
torque_err_initial(iter) = rms(A_B*Phi_B_initial - b);

iter
end
disp('..done')



%% Trajectory Plot
if plot_and_visualize
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
    appVisualizeKUKA(q_optimal);
end