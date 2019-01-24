close all
clear
clc

%%
robot = makeAmbidex();
joints = robot.joints;
motors = robot.motors;

excitation_links = [2,4,6,10];
robot.excitation_links = excitation_links;

partial = false;
if length(excitation_links) ~= joints
    partial = true;
end


A = robot.A;
M = robot.M;

G = zeros(6,6,joints);   % dummy value

num_sample = 10000;

if partial
    cummulativeY = zeros(num_sample*motors, 10*length(excitation_links));
    blocks = zeros(1,length(excitation_links)*10);
    for i = 1:length(excitation_links)
        blocks(i*10-9:i*10) = excitation_links(i)*10-9:excitation_links(i)*10;
    end
else
    cummulativeY = zeros(num_sample*motors, 10*joints + 2*motors);
end


%% Data Generation
for iter = 1:num_sample

%     q = 1000*robot.joint7to10(rand(motors,1));
%     qdot = 1000*robot.joint7to10(rand(motors,1));
%     qddot = 1000*robot.joint7to10(rand(motors,1));

    [q, qdot, qddot] = robot.MotorToJoint(30*rand(motors,1), 10000*rand(motors,1), 10000*rand(motors,1));

    assert(isreal(q));
    assert(isreal(qdot));
    assert(isreal(qddot));
    
    J = robot.motor_jacobian(q);
        
    [tau, V, Vdot] = solveInverseDynamics(A,M,q,qdot,qddot,G, [0;0;0;0;0;9.8]);
    [Y, W] = getRegressorRecursive(A,M,q,V,Vdot);
    
    if partial
        cummulativeY(motors*(iter-1)+1:motors*iter,:) = J'*Y(:,blocks);        
    else
        cummulativeY(motors*(iter-1)+1:motors*iter,:) = [J'*Y getMotorFrictionRegressor(qdot)];
    end
end

[U,S,V] = svd(cummulativeY,'econ');
plot(log(diag(S)))

% B = V(:,1:68 + 2*motors)';
