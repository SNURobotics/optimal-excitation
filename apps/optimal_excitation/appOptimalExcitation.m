%% Trajectory Optimization for Parameter Identification
close all
clear
clc

%% Initialization
robot     = makeWam7();         % robot model
m         = 6;                  % num of trajectory coefs per joints
p_initial = 2*rand(m,robot.dof) - 1;  % initial p

% for i = 1:robot.dof
%     p_initial(:,i) = 0.1*(rand(m,1)-0.5)/robot.Phi(10*(i-1)+1);
% %     p_initial(:,i) = 10/robot.Phi(10*(i-1)+1);
% end

trajectory.order           = 4;      % B Spline cubic base function
trajectory.horizon         = 8;      % trajectory horizon
trajectory.num_sample      = 100;    % number of samples of the trajectory
trajectory.base_frequency  = pi/4;   % Fourier trajectory base frequency

sigma = eye(robot.dof)*1e-5;         % torque covariance
sigma_inv = pinv(sigma);

%% Trajectory Optimization
% w/ metric
options = optimoptions(@fmincon,'Algorithm', 'interior-point', 'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', false, 'TolCon',1e-7,'TolX',1e-5,'MaxFunEvals', ...
    1000000,'MaxIter',1000,'Display','iter','Hessian','bfgs'); %'fin-diff-grads'
% options = optimoptions(@fmincon,'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, 'Algorithm','sqp','TolCon',1e-15,'TolX',1e-15,'MaxFunEvals',1000000,'MaxIter',10000,'Display','iter');
% options = optimoptions(@fmincon, 'Algorithm','sqp','TolCon',1e-15,'TolX',1e-15,'MaxFunEvals',1000000,'MaxIter',10000,'Display','iter');
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp');

tic
[p_optimal,fval,exitflag,output,lam_costate] = fmincon(@(p)getTraceCov(p,robot,trajectory,sigma_inv), p_initial, [], [], [], [], [], [], @(p)getConstraint(p,trajectory), options);
toc

% tic
% [p_optimal,fval,exitflag,output,lam_costate] = fmincon(@(p)getTraceCov(p,robot,trajectory,sigma_inv), p_initial, [], [], [], [], [], [], [], options);
% toc


% p_optimal = p_initial;
% for i = 1 : 39
%     reduced_size = i;
% tic
% [p_optimal,fval,exitflag,output,lam_costate] = fmincon(@(p)getCondNumber(p,robot,trajectory,sigma_inv), p_optimal, [], [], [], [], [], [], [], options);
% toc
% end
% tic
% [p_optimal,fval,exitflag,output,lam_costate] = fmincon(@(p)getrelaxedCondNumber(p,robot,trajectory,sigma_inv), p_initial, [], [], [], [], [], [], [], options);
% toc
% w/o metric
% tic
% [p_baseline,fval_baseline,exitflag_baseline,output_baseline,lam_costate_baseline] = fmincon(@(p)getCondNumberNoMetric(p,robot,trajectory,sigma_inv), p_initial, [], [], [], [], [], [], [], options);
% toc
% 

%% Trajectory Plot
sample_time = linspace(0,trajectory.horizon,trajectory.num_sample);
q_init = makeFourier(p_initial, trajectory.base_frequency, sample_time);
q_opti = makeFourier(p_optimal, trajectory.base_frequency, sample_time);
figure(); hold on;
for i = 1:robot.dof
    subplot(1,robot.dof,i);
    plot(q_init(i,:),'b'); hold on;
    plot(q_opti(i,:),'r');    
end

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