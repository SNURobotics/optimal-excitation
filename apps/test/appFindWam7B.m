close all
clear
clc

%%
robot = wam7robot
n = robot.nDOF;

A = zeros(6,n);
for i=1:n
    A(:,i) = robot.link(i).screw;
end

M = zeros(4,4,n);
for i=1:n
    M(:,:,i) = inverse_SE3(robot.link(i).M);
end

G = zeros(6,6,n);   % dummy value

m = 30;             % number of spline parameters
num_trajectory   = 10;
num_sample       = 100;
horizon          = 30;
trajectory_order = 4;
sample_time      = linspace(0,horizon,num_sample);

cummulativeY = zeros(num_trajectory*num_sample*n, 10*n);

%% Data Generation
for iter = 1:num_trajectory  
    p = rand(m,n)*20*pi;
    % trajectory generation with parameter p
    [q, qdot, qddot]    = makeSpline(p, trajectory_order, horizon, sample_time);

    for t=1:num_sample
        [tau, V, Vdot] = solveInverseDynamics(A,M,q(:,t),qdot(:,t),qddot(:,t),G);
        [Y, W] = getRegressorRecursive(A,M,q(:,t),V,Vdot);
        
        cummulativeY(n*((i-1)*num_sample + (t-1))+1:n*((i-1)*num_sample + t),:) = Y;
    end
end

[U,S,V] = svd(cummulativeY);
B = V(:,1:39)'