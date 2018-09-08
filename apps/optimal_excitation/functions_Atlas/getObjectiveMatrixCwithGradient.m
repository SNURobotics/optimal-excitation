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

    stack = CStack();
    
    % dynamic parameters
    B = robot.B;
    num_base = size(B,1);
    mertic_Phi = robot.pd_metric_Phi;    % metric
    B_metric_inv_Bt = robot.B_metric_inv_Phi_Bt;
    
    G = zeros(6,6,n);   % dummy value
    
    % kinematic variables
    T           = zeros(4,4,n); % T_{i,i-1}
    Ad_T        = zeros(6,6,n); % Ad_T_{i,i-1}
    T_global    = zeros(4,4,n); % T_i0 
    T_global(:,:,robot.root) = eye(4);
    Ad_T_global = zeros(6,6,n); % Ad_T_i0
    Ad_T_global(:,:,robot.root) = eye(6);
    dAd_T_global = zeros(6,6,n,m,n);
    
    
    V_0     = zeros(6,1);       % base velocity
    Vdot_0  = zeros(6,1);       % base acceleration
    dV_0    = zeros(6,1);       % derivative of base velocity
    dVdot_0 = zeros(6,1);       % derivative of base acceleration
    Vdot_0(6) = 9.8;

    V = zeros(6,n);
    Vdot = zeros(6,n);
    dV     = zeros(6,n,m,n);
    dVdot  = zeros(6,n,m,n);
    
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
    dY = zeros(6,10*n,m,n);
    dY_B = zeros(6,num_base,m,n,num_sample);
    sum_A = zeros(num_base, num_base,num_sample);
    sum_dA = zeros(num_base, num_base,m,n,num_sample);
    
    %% Integration
    
    % Spline trajectory generation with parameter p
%     [q, qdot, qddot]    = makeSpline(p, trajectory_order, horizon, sample_time);
%     [dq, dqdot, dqddot] = getSplineDerivative(p, trajectory_order, horizon, sample_time);

    % Spline trajectory generation with parameter p
    [q, qdot, qddot]    = makeFourier(p, base_frequency, sample_time, (robot.q_max + robot.q_min)/2);
    [dq, dqdot, dqddot] = getFourierDerivative(p, base_frequency, sample_time);
    
    for t=1:num_sample
        % get V, Vdot, dV, dVdot from forward recursion      
        stack.push(robot.root);
        while(~stack.isempty)
            i = stack.pop;

            for child = 1:size(robot.tree{i}.children,1)
                stack.push(robot.tree{i}.children(child));
            end

            % T, Ad_T
            T(:,:,i) = exp_se3(-A(:,i)*q(i,t))*M(:,:,i);
            Ad_T(:,:,i) = large_Ad(T(:,:,i));
            T_global(:,:,i) = exp_se3(-robot.A(:,i) * q(i,t))* robot.M(:,:,i) *T_global(:,:,robot.tree{i}.parent);
            Ad_T_global(:,:,i) = large_Ad(T_global(:,:,i));
            
            % dAd_T_global
            for p = 1:m
                for k = 1:n
                    if i ~= robot.root
                        dAd_T_global(:,:,i,p,k) = -dq(i,p,k,t)*small_ad(A(:,i))*Ad_T_global(:,:,i) + Ad_T(:,:,i)*dAd_T_global(:,:,robot.tree{i}.parent,p,k);
                    end
                end
            end

            % V, Vdot
            if i == robot.root
                V(:,i)    = V_0;
                Vdot(:,i) = Vdot_0;
            else
                V(:,i)    = Ad_T(:,:,i)*V(:,robot.tree{i}.parent)    + A(:,i)*qdot(i,t);
                Vdot(:,i) = Ad_T(:,:,i)*Vdot(:,robot.tree{i}.parent) + small_ad(V(:,i))*A(:,i)*qdot(i,t) + A(:,i)*qddot(i,t);
            end
            
            % dV, dVdot
            if i == robot.root
                for p = 1:m
                    for k = 1:n
                        dV(:,i,p,k)    = dV_0;
                        dVdot(:,i,p,k) = dVdot_0;
                    end
                end
            else
                for p = 1:m
                    for k = 1:n
                        dV(:,i,p,k)    = Ad_T(:,:,i)*dV(:,robot.tree{i}.parent,p,k)    - small_ad(A(:,i))*Ad_T(:,:,i)*V(:,robot.tree{i}.parent)*dq(i,p,k,t)    + A(:,i)*dqdot(i,p,k,t);
                        dVdot(:,i,p,k) = Ad_T(:,:,i)*dVdot(:,robot.tree{i}.parent,p,k) - small_ad(A(:,i))*Ad_T(:,:,i)*Vdot(:,robot.tree{i}.parent)*dq(i,p,k,t) + small_ad(dV(:,i,p))*A(:,i)*qdot(i,t) ...
                                        + small_ad(V(:,i))*A(:,i)*dqdot(i,p,k,t) + A(:,i)*dqddot(i,p,k,t);
                    end   
                end
            end            
        end

        for i = 1:n
            % Y
            Y(:,10*(i-1)+1:10*i) = Ad_T_global(:,:,i)'*(convertVelocityToRegressor(Vdot(:,i)) - small_ad(V(:,i))'*convertVelocityToRegressor(V(:,i)));
            
            % dY
            for p=1:m
                for k = 1:n
                    dY(:,10*(i-1)+1:10*i,p,k) = dAd_T_global(:,:,i,p,k)'* (convertVelocityToRegressor(Vdot(:,i)) - small_ad(V(:,i))'*convertVelocityToRegressor(V(:,i))) ...
                                 + Ad_T_global(:,:,i)'*(convertVelocityToRegressor(dVdot(:,i,p,k)) - small_ad(dV(:,i,p,k))'*convertVelocityToRegressor(V(:,i)) ...
                                 - small_ad(V(:,i))'*convertVelocityToRegressor(dV(:,i,p,k)));                    
                end
            end
        end
        
        % Y -> YB
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
%     C = sum(sum_A,3) * B_metric_inv_Bt;
    gradC = zeros(num_base,num_base,m,n); 
    for k = 1:m
        for l = 1:n
%             gradC(:,:,k,l) = sum(sum_dA(:,:,k,l,:), 5) * B_metric_inv_Bt;
            gradC(:,:,k,l) = B_metric_invBt_half* sum(sum_dA(:,:,k,l,:), 5) * B_metric_invBt_half;
            gradC(:,:,k,l) = (gradC(:,:,k,l) + gradC(:,:,k,l)')/2;
        end
    end
end