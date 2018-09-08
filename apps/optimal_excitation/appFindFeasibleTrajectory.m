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

num_trajectory  = 15;  % number of trajectories to find
max_tries = 5000;    % maximum tries

trajectory.order           = 4;       % B Spline cubic base function
trajectory.horizon         = 3;      % trajectory horizon
trajectory.num_sample      = 50;      % number of samples of the trajectory
sample_time      = linspace(0,trajectory.horizon,trajectory.num_sample);

sigma     = eye(6) * 1e-2;            % torque covariance
sigma_inv = pinv(sigma);

m         = 1;                       % num of trajectory coefs per joints
p_initial = rand(m,robot.dof);        % initial p

%% Find Feasible Interpolation Points
qi = zeros(n,num_trajectory+1);
q_total = zeros(n,num_trajectory * trajectory.num_sample);

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
    %         figure(); hold on;
    %         for t = 1:trajectory.num_sample
    %             zmp = getZMP(robot, q(:,t), qdot(:,t), qddot(:,t));
    %             plot(zmp(1), zmp(2), '.'); hold on;
    %         end
    
            q_total(:,1+trajectory.num_sample*(tr-1):trajectory.num_sample*tr) = q;
            p{tr} = p_initial;
            break;
        end
    end
end

%% Visualization
appVisualizeAtlas(q_total)

%% Sample random points within joint limits
function feasible_initial = getInit(robot)
    feasible_initial = zeros(robot.dof,1);
    for i = 1:robot.dof
        feasible_initial(i) = (rand() - 0.5)*(robot.q_max(i) - robot.q_min(i)) + (robot.q_max(i) + robot.q_min(i))/2;
    end
end