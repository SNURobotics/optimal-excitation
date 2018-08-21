close all
clear
clc

%%
robot = makeKukaR820;
n = robot.dof;

A = robot.A;
M = robot.M;

G = zeros(6,6,n);   % dummy value

m = 30;             % number of spline parameters
num_sample       = 1000;

cummulativeY = zeros(num_sample*n, 10*n);

%% Data Generation
for iter = 1:num_sample
    q = 1000*rand(n,1);
    qdot = 1000*rand(n,1);
    qddot = 1000*rand(n,1);

    [tau, V, Vdot] = solveInverseDynamics(A,M,q,qdot,qddot,G);
    [Y, W] = getRegressorRecursive(A,M,q,V,Vdot);
        
    cummulativeY(n*(iter-1)+1:n*iter,:) = Y;
end

[U,S,V] = svd(cummulativeY);
B = V(:,1:39)'