%% Trajectory Optimization for Dynamic Parameter Identification
close all
clear
clc

warning_handler = warning('query','last');
rmpath(genpath('apps/optimal_excitation/functions_hexarotor/'));
rmpath(genpath('apps/optimal_excitation/functions_KUKA/'));
addpath(genpath('apps/optimal_excitation/functions_Atlas/'));
warning('off',warning_handler.identifier);

%% Options
load_trajectory      = true;
visualize_trajectory = false;

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

if ~load_trajectory
    for tr = 1:num_trajectory+1
        while true
            qi(:,tr) = getFeasiblePose(robot);
            com = getCOM(robot, qi(:,tr));
            if com(1) > robot.zmp_x_min && com(1) < robot.zmp_x_max && com(2) > robot.zmp_y_min && com(2) < robot.zmp_y_max 
                break
            end
        end
    end
end

%% Find Feasible Trajectories for the Points
if ~load_trajectory
    for tr = 1:num_trajectory

        tries = 0;
        while true
            tries = tries + 1;

            % if hard to find a feasible trajectory, reset the end point
            if tries > max_tries
                while true
                    qi(:,tr+1) = getFeasiblePose(robot);
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
                 p_initial(i,:) = getFeasiblePose(robot)';
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
end

%% Trajectory Generation from qi & p
if load_trajectory
    load('..\data\qi_atlas_3sec_50samples_100trajectories.mat');
    load('..\data\p_atlas_3sec_50samples_100trajectories.mat');

    for tr = 1:num_trajectory

        [q, qdot, qddot] = makeSplineP2P(qi(:,tr),qi(:,tr+1),p{tr}, trajectory.order, trajectory.horizon, sample_time);

        q_total(:,1+trajectory.num_sample*(tr-1):trajectory.num_sample*tr) = q;
        qdot_total(:,1+trajectory.num_sample*(tr-1):trajectory.num_sample*tr) = qdot;
        qddot_total(:,1+trajectory.num_sample*(tr-1):trajectory.num_sample*tr) = qddot;

    end
end

%% Visualization
if visualize_trajectory
    figure('Name','Atlas v5 :','NumberTitle','off','units','pixels','pos',[-1000 200 900 900]);
    for i = 50
        
        cur = 0;
        for j = 1+trajectory.num_sample*(i-1):8:trajectory.num_sample*i
            cur = cur + 1;
            subplot(2,7,cur);
            appVisualizeAtlas(q_total(:,j));
        end
    end
    for i = 100
        
        for j = 1+trajectory.num_sample*(i-1):7:trajectory.num_sample*i
            cur = cur + 1;
            subplot(2,7,cur);
            appVisualizeAtlas(q_total(:,j));
        end
    end
end

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

pd_metric_Phi_0 = zeros(10*robot.dof, 10*robot.dof);
for i = 1:robot.dof
    pd_metric_Phi_0(10*(i-1)+1:10*i, 10*(i-1)+1:10*i) = getPDMetricInertiaPhi(Phi_0(10*(i-1)+1:10*i));
end

B_metric_inv_Phi_Bt_0 = B * pinv(pd_metric_Phi_0) * B';
B_metric_inv_Phi_Bt_0 = (B_metric_inv_Phi_Bt_0+B_metric_inv_Phi_Bt_0')/2;

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
cum_A = zeros(6*trajectory.num_sample*(num_trajectory/2),num_base);
b = zeros(6*trajectory.num_sample*(num_trajectory/2),1);

% stack for tree traversing
stack = CStack();

disp('computing covariance matrix inverse C..');
% compute C
for t=1:trajectory.num_sample * (num_trajectory/2)
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

B_metric_invBt_half_0 = sqrtm(B_metric_inv_Phi_Bt_0);
B_metric_invBt_half_0 = (B_metric_invBt_half_0+B_metric_invBt_half_0')/2;
% B_metric_invBt_half_0 = eye(size(robot.B_metric_inv_Phi_Bt));

sum_A_sigma_At = (sum_A_sigma_At + sum_A_sigma_At')/2;
C = B_metric_invBt_half_0 * sum_A_sigma_At * B_metric_invBt_half_0;
C = (C+C')/2;

disp('computing Phi_B..');
% eigen decomposition
[Q, D] = eig(C);
D = inv(D);

k = size(D,1);
eigenvalues = zeros(k,1);
for i=1:k
    eigenvalues(i) = D(i,i);
end

V = B_metric_invBt_half_0 * Q;

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

cum_sigma_inv = zeros(6*trajectory.num_sample*(num_trajectory/2));
for t = 1:trajectory.num_sample*(num_trajectory/2)
    cum_sigma_inv(1+ 6*(t-1):6*t,1+ 6*(t-1):6*t) = sigma_inv;
end

Phi_B = V_neg*lambda_neg*V_neg'*cum_A'*cum_sigma_inv*b + V_pos*V_pos'*pinv(B_metric_inv_Phi_Bt_0)*B*Phi_0;
Phi_B_pure = inv(cum_A'*cum_sigma_inv*cum_A)*cum_A'*cum_sigma_inv*b;
disp('done');

temp = [Phi_B B*robot.Phi];

norm((robot.B_metric_inv_Phi_Bt)^(-0.5)*(Phi_B_pure-B*robot.Phi))/sqrt(size(B,1))
figure;
plot(sort(log(eigenvalues)./log(10), 'ascend'),'black'); hold on;
hline =refline([0 -2]);
hline.LineStyle = '--';
hline.Color = 'r';
xlim([1 size(B,1)]);
xlabel('Index of eigenvalue \lambda_i of C_0^{-1/2}CC_{0}^{-1}');
ylabel('log_{10}(\lambda_i)');

disp('done');

%% Prediction Error
disp('Running Prediction Error..')

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

pd_metric_Phi_0 = zeros(10*robot.dof, 10*robot.dof);
for i = 1:robot.dof
    pd_metric_Phi_0(10*(i-1)+1:10*i, 10*(i-1)+1:10*i) = getPDMetricInertiaPhi(Phi_0(10*(i-1)+1:10*i));
end

B_metric_inv_Phi_Bt_0 = B * pinv(pd_metric_Phi_0) * B';
B_metric_inv_Phi_Bt_0 = (B_metric_inv_Phi_Bt_0+B_metric_inv_Phi_Bt_0')/2;

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
cum_A = zeros(6*trajectory.num_sample*(num_trajectory/2),num_base);
b = zeros(6*trajectory.num_sample*(num_trajectory/2),1);

% stack for tree traversing
stack = CStack();

disp('computing covariance matrix inverse C..');
% compute C
for t=trajectory.num_sample * (num_trajectory/2) + 1:trajectory.num_sample * num_trajectory
% for t=1:trajectory.num_sample * (num_trajectory/2)
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
    j = t - (trajectory.num_sample/2) * num_trajectory;
%     j = t;
    cum_A(1+ 6*(j-1):6*j,:) = Y_B;
    
    b(1+ 6*(j-1):6*j,:) = Y*robot.Phi;
end

cum_sigma_inv_half = zeros(6*trajectory.num_sample*(num_trajectory/2));
for t = 1:trajectory.num_sample*(num_trajectory/2)
    cum_sigma_inv_half(1+ 6*(t-1):6*t,1+ 6*(t-1):6*t) = (sigma_inv)^0.5;
end
% [rms(cum_sigma_inv_half*(cum_A * Phi_B - b))]
rms(cum_A([1:6:end, 2:6:end, 3:6:end],:)*B*Phi_0-b([1:6:end, 2:6:end, 3:6:end],:))

disp('done');

%% Joint Torque Estimation
disp('Running Joint Torque Estimation..')

q = q_total; %q = rand(size(q_total))*pi;
qdot = qdot_total; %qdot = rand(size(q_total))*pi;
qddot = qddot_total; %qddot = rand(size(q_total))*pi;

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

% dynamic parameters
B = robot.B;

% stack for tree traversing
stack = CStack();

tau = zeros(robot.dof, trajectory.num_sample * (num_trajectory/2));
tau_predict = zeros(robot.dof, trajectory.num_sample * (num_trajectory/2));
Y_B_tau_all = zeros((robot.dof)*trajectory.num_sample * (num_trajectory/2), size(B,1));


for t=trajectory.num_sample * (num_trajectory/2) + 1:trajectory.num_sample * num_trajectory
    
    t
    
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
    
    Y = getRegressorRecursiveTree(robot, V, Vdot, Ad_T);
    tau(:,t - trajectory.num_sample * (num_trajectory/2)) = Y([1:12,13, 14:end],:)* robot.Phi;
    tau_predict(:,t - trajectory.num_sample * (num_trajectory/2)) = Y([1:12,13, 14:end],:)*B'/(B*B') * Phi_B;
    Y_B_tau_all(1 + (robot.dof)*(t-1-trajectory.num_sample * (num_trajectory/2)) : (robot.dof)*(t-trajectory.num_sample * (num_trajectory/2)), :) = Y([1:12,13, 14:end],:)*B'/(B*B');
end

Phi_B

% rms(tau(:)-tau_predict(:))
% rms(Y_B_tau_all*Phi_B - tau(:))

% rms(Y_B_tau_all*(Phi_B-B*robot.Phi))
figure(1); hold on;
plot((Y_B_tau_all*(Phi_B - B*robot.Phi)));
figure(2);
plot((Y_B_tau_all*(Phi_B_pure - B*robot.Phi)));


% rms(sqrt(((Y_B_tau_all*(Phi_B_pure - B*robot.Phi)).^2)./((Y_B_tau_all*(Phi_B - B*robot.Phi)).^2)))
err_ratio = sqrt(((Y_B_tau_all*(Phi_B_pure - B*robot.Phi)).^2)./((Y_B_tau_all*(Phi_B - B*robot.Phi)).^2));
figure(1);cur = 0;
for i =[1:12, 14: robot.dof]
    cur = cur + 1;
    err_dof(i) = rms(Y_B_tau_all([i:(robot.dof):end],:)*(Phi_B - B*robot.Phi))/rms(Y_B_tau_all([i:(robot.dof):end],:)*(B*robot.Phi));
    err_dof_pure(i) = rms(Y_B_tau_all([i:(robot.dof):end],:)*(Phi_B_pure - B*robot.Phi))/rms(Y_B_tau_all([i:(robot.dof):end],:)*(B*robot.Phi));
%     err_ratio_dof(i) = rms(err_ratio(i : (robot.dof-1):end,1));
    subplot(3,10,cur); hold on;
    plot(linspace(0, 150, size(Y_B_tau_all,1)/(robot.dof)), Y_B_tau_all([i:(robot.dof):end],:)*(Phi_B_pure - B*robot.Phi), 'r');
    plot(linspace(0, 150, size(Y_B_tau_all,1)/(robot.dof)), Y_B_tau_all([i:(robot.dof):end],:)*(Phi_B - B*robot.Phi), 'b');
    
    hline =refline([0 0]);
    hline.LineStyle = '--';
    hline.Color = 'black';
    xlim([0 150]);
end
legend('w/o param. selection', 'w/ param. selection');

figure;
bar([err_dof;err_dof_pure]');
temp = [err_dof./err_dof_pure];
temp = [err_dof;err_dof_pure]'*100;

disp('done');
