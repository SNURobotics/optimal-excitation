%% Trajectory Optimization for Dynamic Parameter Identification
close all
clear
clc

rmpath(genpath('apps/optimal_excitation/functions_hexarotor/'));
rmpath(genpath('apps/optimal_excitation/functions/'));
addpath(genpath('apps/optimal_excitation/functions_Atlas/'));

%% Initialization
disp('initializing..')

robot     = makeAtlasV5();           % robot model

n = robot.dof;

num_trajectory  = 100;  % number of trajectories to find
max_tries = 5000;    % maximum tries

trajectory.order           = 4;       % B Spline cubic base function
trajectory.horizon         = 3;       % trajectory horizon
trajectory.num_sample      = 50;      % number of samples of the trajectory
sample_time      = linspace(0,trajectory.horizon,trajectory.num_sample);

sigma     = diag([0.25 0.25 0.25 0.5 0.5 0.5].^2);            % force covariance
sigma_inv = pinv(sigma);

m         = 1;                        % num of trajectory coefs per joints
p_initial = rand(m,robot.dof);        % initial p

%% Find Feasible Interpolation Points
qi = zeros(n,num_trajectory+1);
q_total = zeros(n,num_trajectory * trajectory.num_sample);
qdot_total = zeros(n,num_trajectory * trajectory.num_sample);
qddot_total = zeros(n,num_trajectory * trajectory.num_sample);

for tr = 1:num_trajectory+1
    while true
        qi(:,tr) = getInit(robot);
        com = getCOM(robot, qi(:,tr));
        if com(1) > robot.zmp_x_min && com(1) < robot.zmp_x_max && com(2) > robot.zmp_y_min && com(2) < robot.zmp_y_max 
            break
        end
    end
end

%% Find Feasible Trajectories for the Points
for tr = 1:num_trajectory

    tries = 0;
    while true
        tries = tries + 1;
        
        % if hard to find a feasible trajectory, reset the end point
        if tries > max_tries
            while true
                qi(:,tr+1) = getInit(robot);
                com = getCOM(robot, qi(:,tr+1));
                if com(1) > robot.zmp_x_min && com(1) < robot.zmp_x_max && com(2) > robot.zmp_y_min && com(2) < robot.zmp_y_max 
                    break
                end
            end
            tries = 1;
        end
        
        disp([num2str(tr) 'th trajectory ' num2str(tries) 'th try']);

         p_initial = zeros(m,robot.dof);
         for i = 1:m
             p_initial(i,:) = getInit(robot)';
         end
        [q, qdot, qddot] = makeSplineP2P(qi(:,tr),qi(:,tr+1),p_initial, trajectory.order, trajectory.horizon, sample_time);

        flag = true;
        for t = 1:trajectory.num_sample
            if size(find(qdot(:,t) - robot.qdot_min > 0),1) ~= n || size(find(qdot(:,t) - robot.qdot_max < 0),1) ~= n
                flag = false;
                break;
            end    
            zmp = getZMP(robot, q(:,t), qdot(:,t), qddot(:,t));
            if zmp(1) < robot.zmp_x_min || zmp(1) > robot.zmp_x_max || zmp(2) < robot.zmp_y_min || zmp(2) > robot.zmp_y_max
                flag = false;
                break;
            end
        end

        if flag == true
            % plot CoM trajectory
    %         figure(); hold on;
    %         for t = 1:trajectory.num_sample
    %             zmp = getZMP(robot, q(:,t), qdot(:,t), qddot(:,t));
    %             plot(zmp(1), zmp(2), '.'); hold on;
    %         end
    
            q_total(:,1+trajectory.num_sample*(tr-1):trajectory.num_sample*tr) = q;
            qdot_total(:,1+trajectory.num_sample*(tr-1):trajectory.num_sample*tr) = qdot;
            qddot_total(:,1+trajectory.num_sample*(tr-1):trajectory.num_sample*tr) = qddot;
            p{tr} = p_initial;
            break;
        end
    end
end

%% Trajectory Generation from qi & p
load('..\data\qi_atlas_3sec_50samples_100trajectories.mat');
load('..\data\p_atlas_3sec_50samples_100trajectories.mat');

for tr = 1:num_trajectory

    [q, qdot, qddot] = makeSplineP2P(qi(:,tr),qi(:,tr+1),p{tr}, trajectory.order, trajectory.horizon, sample_time);

    q_total(:,1+trajectory.num_sample*(tr-1):trajectory.num_sample*tr) = q;
    qdot_total(:,1+trajectory.num_sample*(tr-1):trajectory.num_sample*tr) = qdot;
    qddot_total(:,1+trajectory.num_sample*(tr-1):trajectory.num_sample*tr) = qddot;

end

%% Visualization
appVisualizeAtlas(q_total)


%% Parameter Identification
disp('Running Parameter Identification..')

q = q_total;
qdot = qdot_total;
qddot = qddot_total;

% dynamic parameters
B = robot.B;
num_base = size(B,1);

Phi_0 = zeros(10*robot.dof,1);  % nominal phi
for i = 1:robot.dof
    Phi_0(1+(i-1)*10:i*10) = getNominalPhi(robot.Phi(1+(i-1)*10:i*10));
end

pd_metric_Phi = zeros(10*robot.dof, 10*robot.dof);
for i = 1:robot.dof
    pd_metric_Phi(10*(i-1)+1:10*i, 10*(i-1)+1:10*i) = getPDMetricInertiaPhi(Phi_0(10*(i-1)+1:10*i));
end

B_metric_inv_Phi_Bt = B * pinv(pd_metric_Phi) * B';
B_metric_inv_Phi_Bt = (B_metric_inv_Phi_Bt+B_metric_inv_Phi_Bt')/2;

% kinematic parameters
A = robot.A;
M = robot.M;

% kinematic variables
T           = zeros(4,4,n);         % T_{i,i-1}
Ad_T        = zeros(6,6,n);         % Ad_T_{i,i-1}
T_global    = zeros(4,4,n);         % T_i0 
T_global(:,:,robot.root) = eye(4);
Ad_T_global = zeros(6,6,n);         % Ad_T_i0
Ad_T_global(:,:,robot.root) = eye(6);

V_0     = zeros(6,1);               % base velocity
Vdot_0  = zeros(6,1);               % base acceleration
Vdot_0(6) = 9.8;

V = zeros(6,n);
Vdot = zeros(6,n);

% size initialization
Y = zeros(6, 10*robot.dof);
sum_A_sigma_At = zeros(num_base, num_base);
cum_A = zeros(6*trajectory.num_sample*num_trajectory,num_base);
b = zeros(6*trajectory.num_sample*num_trajectory,1);

% stack for tree traversing
stack = CStack();

% compute C
for t=1:trajectory.num_sample * num_trajectory
    % get V, Vdot, dV, dVdot from forward recursion      
    stack.push(robot.root);
    while(~stack.isempty)
        i = stack.pop;

        for child = 1:size(robot.tree{i}.children,1)
            stack.push(robot.tree{i}.children(child));
        end

        % T, Ad_T
        T(:,:,i) = exp_se3(-A(:,i)*q(i,t))*M(:,:,i);
        Ad_T(:,:,i) = large_Ad(T(:,:,i));
        T_global(:,:,i) = exp_se3(-A(:,i) * q(i,t))* M(:,:,i) *T_global(:,:,robot.tree{i}.parent);
        Ad_T_global(:,:,i) = large_Ad(T_global(:,:,i));

        % V, Vdot
        if i == robot.root
            V(:,i)    = V_0;
            Vdot(:,i) = Vdot_0;
        else
            V(:,i)    = Ad_T(:,:,i)*V(:,robot.tree{i}.parent)    + A(:,i)*qdot(i,t);
            Vdot(:,i) = Ad_T(:,:,i)*Vdot(:,robot.tree{i}.parent) + small_ad(V(:,i))*A(:,i)*qdot(i,t) + A(:,i)*qddot(i,t);
        end
    end

    for i = 1:n
        % Y
        Y(:,10*(i-1)+1:10*i) = Ad_T_global(:,:,i)'*(convertVelocityToRegressor(Vdot(:,i)) - small_ad(V(:,i))'*convertVelocityToRegressor(V(:,i)));
    end

    % Y -> YB
    Y_B = Y*B'/(B*B');
    sum_A_sigma_At = sum_A_sigma_At + Y_B' * sigma_inv * Y_B;
    cum_A(1+ 6*(t-1):6*t,:) = Y_B;
    
    b(1+ 6*(t-1):6*t,:) = Y*robot.Phi + sqrtm(sigma)*randn(6,1);
end

B_metric_invBt_half = sqrtm(B_metric_inv_Phi_Bt);
B_metric_invBt_half = (B_metric_invBt_half+B_metric_invBt_half')/2;

sum_A_sigma_At = (sum_A_sigma_At + sum_A_sigma_At')/2;
C = B_metric_invBt_half * sum_A_sigma_At * B_metric_invBt_half;
C = (C+C')/2;
disp('covariance matrix inverse C has been computed');

% eigen decomposition
[Q, D] = eig(C);
D = inv(D);

k = size(D,1);
eigenvalues = zeros(k,1);
for i=1:k
    eigenvalues(i) = D(i,i);
end

V = B_metric_invBt_half * Q;

max_variance = 1e-2;
num_identify = size(find(eigenvalues < max_variance),1);

V_neg = zeros(size(V,1), num_identify);
V_pos = zeros(size(V,1), size(V,1) - num_identify);
lambda_neg = zeros(num_identify); 
lambda_pos = zeros(size(V,1) - num_identify);

neg = 0;
pos = 0;
for i = 1:size(V,1)
    if D(i,i) < max_variance
        neg = neg + 1;
        V_neg(:,neg) = V(:,i);
        lambda_neg(neg,neg) = D(i,i);
    else
        pos = pos + 1;
        V_pos(:,pos) = V(:,i);
        lambda_pos(pos,pos) = D(i,i);
    end
end

cum_sigma_inv = zeros(6*trajectory.num_sample*num_trajectory);
for t = 1:trajectory.num_sample*num_trajectory
    cum_sigma_inv(1+ 6*(t-1):6*t,1+ 6*(t-1):6*t) = sigma_inv;
end

Phi_B = V_neg*lambda_neg*V_neg'*cum_A'*cum_sigma_inv*b + V_pos*V_pos'*pinv(B_metric_inv_Phi_Bt)*B*Phi_0;

disp('done');
[pinv(B_metric_invBt_half)*Phi_B pinv(B_metric_invBt_half)*B*robot.Phi];