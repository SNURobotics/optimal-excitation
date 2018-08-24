%% Compute the Objective Matrix C = sum[Y_B'*Sigma^-1*Y_B]B*G(phi_0)*B'
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
    
    G = zeros(6,6,n);                    % dummy value
    Vdot_0 = [0; 0; 0; 0; 0; 9.8];      % gravity

    % kinematic parameters
    A = robot.A;
    M = robot.M;

    S_inv = robot.S_inv;

    % parameterize q, qdot, qddot itself
    num_sample = 1;
    
    % Spline trajectory parameters
%     num_sample       = trajectory.num_sample;
%     horizon          = trajectory.horizon;
%     trajectory_order = trajectory.order;
%     sample_time      = linspace(0,horizon,num_sample);
%     sample_interval  = horizon/num_sample;              % dt for numerical integration

    % Fourier trajectory parameters
%     num_sample       = trajectory.num_sample;
%     horizon          = trajectory.horizon;
%     base_frequency   = trajectory.base_frequency;
%     sample_time      = linspace(0,horizon,num_sample);
%     sample_interval  = horizon/num_sample;                % dt for numerical integration

    % size initialization
    sum_A = zeros(num_base, num_base,num_sample);

    %% Integration
    
    % Spline trajectory generation with parameter p
%     [q, qdot, qddot]    = makeSpline(p, trajectory_order, horizon, sample_time);

    % Fourier trajectory generation with parameter p
%     [q, qdot, qddot]    = makeFourier(p, base_frequency, sample_time);
    
    % parameterize q, qdot, qddot itself
    q = p(1,:)';
    qdot = p(2,:)';
    qddot = p(3,:)';  
    dq = zeros(n,m,n);
    dqdot= zeros(n,m,n);
    dqddot = zeros(n,m,n);
    for i = 1:n
        dq(i,1,i) = 1;
        dqdot(i,2,i) = 1;
        dqddot(i,3,i) = 1;
    end
    
    for t=1:num_sample
        % get V, Vdot, dV, dVdot from forward recursion, use dummy G and F      
        [tau, V, Vdot, F] = solveInverseDynamics(A,M,q(:,t),qdot(:,t),qddot(:,t),G,Vdot_0);

        % get Y, dY by recursion
        Y = S_inv*(convertVelocityToRegressor(Vdot(:,n)) - small_ad(V(:,n))'*convertVelocityToRegressor(V(:,n)));
        Y_B = Y*B'*pinv(B*B');
        sum_A(:,:,t) = Y_B'*sigma_inv*Y_B;
    end
    
    C = sum(sum_A,3) * B_metric_inv_Bt;
end