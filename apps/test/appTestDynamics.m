close all
clear
clc

%%
robot = makeKukaR820;

n = robot.dof;
m = 3;
A = robot.A;
M = robot.M;
G = robot.G;
Phi = robot.Phi;

q = rand(7,1);
qdot = rand(7,1);
qddot = rand(7,1);

%% Inverse Dynamics
[tau_from_inverse_dynamics, V, Vdot, F] = solveInverseDynamics(A,M,q,qdot,qddot,G);

[Y, W] = getRegressorRecursive(A,M,q,V,Vdot);

tau_from_regressor        = Y*Phi;
tau = [tau_from_inverse_dynamics, tau_from_regressor]

dq = rand(n,m,n);
dqdot = rand(n,m,n);
dqddot = rand(n,m,n);

[dtau, dV, dVdot] = solveInverseDynamicsDerivatives(A,M,q,qdot,G,V,Vdot,dq,dqdot,dqddot,F);
dY = getRegressorDerivativesRecursive(A,M,q,V,Vdot,dq,dV,dVdot,W);

%% Forward Dynamics
qddot_f = solveForwardDynamics(A,M,q,qdot,tau,G)