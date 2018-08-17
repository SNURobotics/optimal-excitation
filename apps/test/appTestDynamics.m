close all
clear
clc

%%
A = [0, 0, 1, 0, 1, 0; 1, 0, 0, 1, 0, 0; 1,1,0,0,0,0]';
M = zeros(4,4,3);
M(1:3,1:3,1) = eye(3);
M(1:3,1:3,2) = eye(3);
M(1:3,1:3,3) = eye(3);
M(1:3,4,1) = [1, 4, 5]';
M(1:3,4,2) = [1, 0, 0]';
M(1:3,4,3) = [1, -1,6]';
M(4,4,:) = 1;

q = [1, 3, -2]';
qdot = [3, 4, 1]';
qddot = [3, 5, 8]';

G(:,:,1) = 6*eye(6);
G(:,:,2) = 1*eye(6);
G(:,:,3) = 1*eye(6);

%% Inverse Dynamics
[tau, V, Vdot, F] = solveInverseDynamics(A,M,q,qdot,qddot,G);

[Y, W] = getRegressorRecursive(A,M,q,V,Vdot);
phi(1:10,1) = convertInertiaGToPhi(G(:,:,1));
phi(11:20,1) = convertInertiaGToPhi(G(:,:,2));
phi(21:30,1) = convertInertiaGToPhi(G(:,:,3));

tau_from_inverse_dynamics = tau
tau_from_regressor        = Y*phi

dq = rand(3,4,3);
dqdot = dq;
dqddot= dq;

[dtau, dV, dVdot] = solveInverseDynamicsDerivatives(A,M,q,qdot,G,V,Vdot,dq,dqdot,dqddot,F);
dY = getRegressorDerivativesRecursive(A,M,q,V,Vdot,dq,dV,dVdot,W);

%% Forward Dynamics
qddot_f = solveForwardDynamics(A,M,q,qdot,tau,G)