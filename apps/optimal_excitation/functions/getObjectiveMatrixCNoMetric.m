%% Compute the Objective Matrix C = sum[Y_B'*Sigma^-1*Y_B] and Its Derivative dC/dP
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]      [Description]                                                [Size]
%  p           trajectory parameter to optimize                             m*n
%  robot       robot object (A, B, M)                                       struct
%  trajectory  object w/ trajectory parameters(order, horizon, num_sample)  struct
%  sigma_inv   inverse of torque covariance matrix                          n*n

%% Outputs

%% Implementation
function [C, gradC] = getObjectiveMatrixCNoMetric(p, robot, trajectory, sigma_inv)
    %% Initialization
    n = robot.dof;        % number of joints
    m = size(p,1);        % number of parameters

    % dynamic parameters
    B = robot.B;
    num_base = size(B,1);
    
    G = zeros(6,6,n);   % dummy value
    
    % kinematic parameters
    A = robot.A;
    M = robot.M;
        
    % trajectory parameters
    num_sample       = trajectory.num_sample;
    horizon          = trajectory.horizon;
    trajectory_order = trajectory.order;
    sample_time      = linspace(0,horizon,num_sample);
    sample_interval  = horizon/num_sample;              % dt for numerical integration

    % size initialization
    dY_B = zeros(n,num_base,m);
    sum_A = zeros(num_base, num_base);
    sum_dA = zeros(num_base, num_base,m);
    
    %% Integration
    
    % trajectory generation with parameter p
    [q, qdot, qddot]    = makeSpline(p, trajectory_order, horizon, sample_time);
    [dq, dqdot, dqddot] = getSplineDerivative(p, trajectory_order, horizon, sample_time);

    for t=1:num_sample
        % get V, Vdot, dV, dVdot from forward recursion, use dummy G and F
        [tau, V, Vdot, F] = solveInverseDynamics(A,M,q(:,t),qdot(:,t),qddot(:,t),G);
        [dtau, dV, dVdot] = solveInverseDynamicsDerivatives(A,M,q(:,t),qdot(:,t),G,V,Vdot,dq(:,:,t),dqdot(:,:,t),dqddot(:,:,t),F);

        % get Y, dY by recursion
        [Y, W] = getRegressorRecursive(A,M,q(:,t),V,Vdot);
        dY = getRegressorDerivativesRecursive(A,M,q(:,t),V,Vdot,dq(:,:,t),dV,dVdot,W);
        
        Y_B = Y*B'*pinv(B*B');
        sum_A = sum_A + Y_B'*sigma_inv*Y_B;

        for k = 1:m
            dY_B(:,:,k) = dY(:,:,k)*B'*pinv(B*B');
            sum_dA(:,:,k) = sum_dA(:,:,k) + dY_B(:,:,k)'*sigma_inv*Y_B + Y_B'*sigma_inv*dY_B(:,:,k);
        end
    end
    
    C = sum_A * sample_interval;
    gradC = sum_dA * sample_interval;
end