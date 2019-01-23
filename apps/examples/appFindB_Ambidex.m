close all
clear
clc

%%
robot = makeAmbidex();
joints = robot.joints;
motors = robot.motors;

A = robot.A;
M = robot.M;

G = zeros(6,6,joints);   % dummy value

num_sample = 1000;

cummulativeY = zeros(num_sample*motors, 10*joints);

%% Data Generation
for iter = 1:num_sample
    q = 1000*robot.joint7to10(rand(motors,1));
    qdot = 1000*robot.joint7to10(rand(motors,1));
    qddot = 1000*robot.joint7to10(rand(motors,1));

    J = robot.motor_jacobian(q);
        
    [tau, V, Vdot] = solveInverseDynamics(A,M,q,qdot,qddot,G, [0;0;0;0;0;9.8]);
    [Y, W] = getRegressorRecursive(A,M,q,V,Vdot);
        
    cummulativeY(motors*(iter-1)+1:motors*iter,:) = J'*Y;
end

[U,S,V] = svd(cummulativeY,'econ');
B = V(:,1:68)'