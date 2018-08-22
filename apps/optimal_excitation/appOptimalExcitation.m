%% Trajectory Optimization for Dynamic Parameter Identification
close all
clear
clc

%% Initialization
robot     = makeKukaR820();           % robot model
m         = 10;                       % num of trajectory coefs per joints
p_initial = 2*rand(m,robot.dof) - 1;  % initial p

% for i = 1:robot.dof
%     p_initial(:,i) = 0.1*(rand(m,1)-0.5)/robot.Phi(10*(i-1)+1);
% %     p_initial(:,i) = 10/robot.Phi(10*(i-1)+1);
% end

trajectory.order           = 4;       % B Spline cubic base function
trajectory.horizon         = 20;      % trajectory horizon
trajectory.num_sample      = 10;      % number of samples of the trajectory
trajectory.base_frequency  = pi*0.1;  % Fourier trajectory base frequency

sigma     = eye(robot.dof);           % torque covariance
sigma_inv = pinv(sigma);

%% Trajectory Optimization
% w/ metric
options = optimoptions(@fmincon,'Algorithm', 'interior-point', 'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', false, 'TolCon',1e-7,'TolX',1e-5,'MaxFunEvals', ...
    1000000,'MaxIter',1000,'Display','iter','Hessian','bfgs'); %'fin-diff-grads'
% options = optimoptions(@fmincon,'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, 'Algorithm','sqp','TolCon',1e-15,'TolX',1e-15,'MaxFunEvals',1000000,'MaxIter',10000,'Display','iter');
% options = optimoptions(@fmincon, 'Algorithm','sqp','TolCon',1e-15,'TolX',1e-15,'MaxFunEvals',1000000,'MaxIter',10000,'Display','iter');
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp');

tic
[p_optimal,fval,exitflag,output,lam_costate] = fmincon(@(p)getTraceCov(p,robot,trajectory,sigma_inv), p_initial, [], [], [], [], [], [], @(p)getConstraint(p,trajectory,robot), options);
toc

% tic
% [p_optimal,fval,exitflag,output,lam_costate] = fmincon(@(p)getTraceCov(p,robot,trajectory,sigma_inv), p_initial, [], [], [], [], [], [], [], options);
% toc

%% Least-Square Parameter Identification
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