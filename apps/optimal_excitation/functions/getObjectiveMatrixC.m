%% Compute the Objective Matrix C = sum[Y_B'*Sigma^-1*Y_B]B*G(phi_0)*B' and Its Derivative dC/dP
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]      [Description]                                                [Size]
%  p           trajectory parameter to optimize                             m*n
%  robot       robot object (A, B, M)                                       struct
%  trajectory  object w/ trajectory parameters(order, horizon, num_sample)  struct
%  sigma_inv   inverse of torque covariance matrix                          n*n

%% Outputs
% [Name]      [Description]               [Size]
%  C           objective matrix C          num_base*num_base                                          

%% Implementation
function C = getObjectiveMatrixC(p, robot, trajectory, sigma_inv)
    %% Initialization
    n = robot.dof;        % number of joints
    m = size(p,1);        % number of parameters

    % dynamic parameters
    B = robot.B;
    num_base = size(B,1);
    mertic_Phi = robot.pd_metric_Phi;    % metric
    B_metric_inv_Bt = robot.B_metric_inv_Phi_Bt;
    
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
    sum_A = zeros(num_base, num_base,num_sample);

    %% Integration
    
    % trajectory generation with parameter p
    [q, qdot, qddot]    = makeSpline(p, trajectory_order, horizon, sample_time);
    [dq, dqdot, dqddot] = getSplineDerivative(p, trajectory_order, horizon, sample_time);

    parfor t=1:num_sample
        % get V, Vdot, dV, dVdot from forward recursion, use dummy G and F      
        [tau, V, Vdot, F] = solveInverseDynamics(A,M,q(:,t),qdot(:,t),qddot(:,t),G);

        % get Y, dY by recursion
        [Y, W] = getRegressorRecursive(A,M,q(:,t),V,Vdot);
        Y_B = Y*B'*pinv(B*B');
        sum_A(:,:,t) = Y_B'*sigma_inv*Y_B;
    end
    
    C = sum(sum_A,3) * B_metric_inv_Bt * sample_interval;
end