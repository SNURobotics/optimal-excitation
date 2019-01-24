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

%% Implementation
function C = getObjectiveMatrixC(p, robot, trajectory, sigma_inv)
    %% Initialization
    motors = robot.motors;  % number of motors
    joints = robot.joints;  % number of joints
    m = size(p,1);          % number of parameters per joint

    excitation_links = robot.excitation_links;

    partial = false;
    if length(excitation_links) ~= joints
        partial = true;
        
        blocks = zeros(1,length(excitation_links)*10);
        for i = 1:length(excitation_links)
            blocks(i*10-9:i*10) = excitation_links(i)*10-9:excitation_links(i)*10;
        end
    end
    
    % dynamic parameters
    B = robot.B;
    num_base = size(B,1);
    B_metric_inv_Bt = robot.B_metric_inv_Phi_Bt_0;
    
    G = zeros(6,6,joints);   % dummy value
    
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
    sample_interval  = horizon/num_sample;              % dt for numerical integration

    % size initialization
    sum_A = zeros(num_base, num_base,num_sample);

    J = zeros(joints, motors, num_sample);

    %% Integration
    
     % Spline trajectory generation with parameter p
%     [q, qdot, qddot]    = makeSpline(p, trajectory_order, horizon, sample_time);

    % Spline trajectory generation with parameter p
    [q, qdot, qddot]    = makeFourier(p, base_frequency, sample_time);

    q = robot.joint7to10(q);
    qdot = robot.joint7to10(qdot);
    qddot = robot.joint7to10(qddot);
    
    for t=1:num_sample
        J(:,:,t) = robot.motor_jacobian(q(:,t));
    end
    
    parfor t=1:num_sample
        % get V, Vdot, dV, dVdot from forward recursion, use dummy G and F      
        [tau, V, Vdot, F] = solveInverseDynamics(A,M,q(:,t),qdot(:,t),qddot(:,t),G, [0;0;0;0;0;9.8]);

        % get Y, dY by recursion
        [Y, W] = getRegressorRecursive(A,M,q(:,t),V,Vdot);
         if partial
            Y_B = J(:,:,t)'*Y(:,blocks)*B'*pinv(B*B');
        else
            Y_B = J(:,:,t)'*Y*B'*pinv(B*B');
        end        
        sum_A(:,:,t) = Y_B'*sigma_inv*Y_B;
    end
    
    B_metric_invBt_half = sqrtm(B_metric_inv_Bt);
    B_metric_invBt_half = (B_metric_invBt_half+B_metric_invBt_half')/2;

    C = B_metric_invBt_half* sum(sum_A,3) * B_metric_invBt_half;
    C = (C+C')/2;
    
end