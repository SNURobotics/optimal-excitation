%% Compute the Objective Matrix C = sum[Y_B'*Sigma^-1*Y_B]B*G(phi_0)*B' and Its Derivative dC/dP
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]      [Description]                                                [Size]
%  p           trajectory parameter to optimize                             m*n
%  robot       robot object (A, B, M)                                       struct
%  trajectory  object w/ trajectory parameters                              struct
%  sigma_inv   inverse of torque covariance matrix                          n*n

%% Outputs
% [Name]      [Description]               [Size]
%  C           objective matrix C          num_base*num_base                                          
%  gradC       gradient of C               num_base*num_base*m*n

%% Implementation
function [C, gradC] = getObjectiveMatrixCwithGradient(p, robot, trajectory, sigma_inv)
    %% Initialization
    n = robot.dof;        % number of joints
    m = size(p,1);        % number of parameters

    % dynamic parameters
    B = robot.B;
    num_base = size(B,1);
    B_metric_inv_Bt = robot.B_metric_inv_Phi_Bt_0;
    
    G = zeros(6,6,n);   % dummy value
    
    % kinematic parameters
    A = robot.A;
    M = robot.M;
        
    % Spline trajectory parameters
%     num_sample       = trajectory.num_sample;
%     horizon          = trajectory.horizon;
%     trajectory_order = trajectory.order;
%     sample_time      = linspace(0,horizon,num_sample);
%     sample_interval  = horizon/num_sample;              % dt for numerical integration

    % Fourier trajectory parameters
    num_sample       = trajectory.num_sample;
    horizon          = trajectory.horizon;
    base_frequency   = trajectory.base_frequency;
    sample_time      = linspace(0,horizon,num_sample);
    sample_interval  = horizon/num_sample;                % dt for numerical integration

    % size initialization
    dY_B = zeros(n,num_base,m,n,num_sample);
    sum_A = zeros(num_base, num_base,num_sample);
    sum_dA = zeros(num_base, num_base,m,n,num_sample);
    
    %% Integration
    
    % Spline trajectory generation with parameter p
%     [q, qdot, qddot]    = makeSpline(p, trajectory_order, horizon, sample_time);
%     [dq, dqdot, dqddot] = getSplineDerivative(p, trajectory_order, horizon, sample_time);

    % Fourier trajectory generation with parameter p
    [q, qdot, qddot]    = makeFourier(p, base_frequency, sample_time);
    [dq, dqdot, dqddot] = getFourierDerivative(p, base_frequency, sample_time);

    parfor t=1:num_sample
        % get V, Vdot, dV, dVdot from forward recursion, use dummy G and F      
        [tau, V, Vdot, F] = solveInverseDynamics(A,M,q(:,t),qdot(:,t),qddot(:,t),G);
        [dtau, dV, dVdot] = solveInverseDynamicsDerivatives(A,M,q(:,t),qdot(:,t),G,V,Vdot,dq(:,:,:,t),dqdot(:,:,:,t),dqddot(:,:,:,t),F);

        % get Y, dY by recursion
        [Y, W] = getRegressorRecursive(A,M,q(:,t),V,Vdot);
        dY = getRegressorDerivativesRecursive(A,M,q(:,t),V,Vdot,dq(:,:,:,t),dV,dVdot,W);
        Y_B = Y*B'*pinv(B*B');
        sum_A(:,:,t) = Y_B'*sigma_inv*Y_B;

        for k = 1:m
            for l = 1:n
                dY_B(:,:,k,l,t) = dY(:,:,k,l)*B'*pinv(B*B');
                sum_dA(:,:,k,l,t) = dY_B(:,:,k,l,t)'*sigma_inv*Y_B + Y_B'*sigma_inv*dY_B(:,:,k,l,t);
            end
        end
    end
    
    B_metric_invBt_half = sqrtm(B_metric_inv_Bt);
    B_metric_invBt_half = (B_metric_invBt_half+B_metric_invBt_half')/2;

    C = B_metric_invBt_half* sum(sum_A,3) * B_metric_invBt_half;
    C = (C+C')/2;
    
    gradC = zeros(num_base,num_base,m,n); 
    for k = 1:m
        for l = 1:n
            gradC(:,:,k,l) = B_metric_invBt_half* sum(sum_dA(:,:,k,l,:), 5) * B_metric_invBt_half;
        end
    end
end