%% Compute the Objective Matrix C of P2P Spline Trajectory
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]      [Description]                                                [Size]
%  qi          initial pose                                                 n*1
%  qf          final pose                                                   n*1
%  p           trajectory parameter to optimize                             m*n
%  robot       robot object (A, B, M)                                       struct
%  trajectory  object w/ trajectory parameters                              struct
%  sigma_inv   inverse of torque covariance matrix                          n*n

%% Outputs
% [Name]      [Description]               [Size]
%  C           objective matrix C          num_base*num_base                                          

%% Implementation
function C = getObjectiveMatrixC_Spline(qi, qf, p, robot, trajectory, sigma_inv)
    n = robot.dof;
    sample_time      = linspace(0,trajectory.horizon,trajectory.num_sample);

    [q, qdot, qddot] = makeSplineP2P(qi,qf,p, trajectory.order, trajectory.horizon, sample_time);

    % dynamic parameters
    B = robot.B;
    num_base = size(B,1);
    B_metric_inv_Bt = robot.B_metric_inv_Phi_Bt;

    % kinematic parameters
    A = robot.A;
    M = robot.M;

    % kinematic variables
    T           = zeros(4,4,n);         % T_{i,i-1}
    Ad_T        = zeros(6,6,n);         % Ad_T_{i,i-1}
    T_global    = zeros(4,4,n);         % T_i0 
    T_global(:,:,robot.root) = eye(4);
    Ad_T_global = zeros(6,6,n);         % Ad_T_i0
    Ad_T_global(:,:,robot.root) = eye(6);

    V_0     = zeros(6,1);               % base velocity
    Vdot_0  = zeros(6,1);               % base acceleration
    Vdot_0(6) = 9.8;

    V = zeros(6,n);
    Vdot = zeros(6,n);

    % size initialization
    Y = zeros(6, 10*robot.dof);
    sum_A = zeros(num_base, num_base);

    % stack for tree traversing
    stack = CStack();

    % compute C
    for t=1:trajectory.num_sample
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
            T_global(:,:,i) = exp_se3(-A(:,i) * q(i,t))* M(:,:,i) *T_global(:,:,robot.tree{i}.parent);
            Ad_T_global(:,:,i) = large_Ad(T_global(:,:,i));

            % V, Vdot
            if i == robot.root
                V(:,i)    = V_0;
                Vdot(:,i) = Vdot_0;
            else
                V(:,i)    = Ad_T(:,:,i)*V(:,robot.tree{i}.parent)    + A(:,i)*qdot(i,t);
                Vdot(:,i) = Ad_T(:,:,i)*Vdot(:,robot.tree{i}.parent) + small_ad(V(:,i))*A(:,i)*qdot(i,t) + A(:,i)*qddot(i,t);
            end
        end

        for i = 1:n
            % Y
            Y(:,10*(i-1)+1:10*i) = Ad_T_global(:,:,i)'*(convertVelocityToRegressor(Vdot(:,i)) - small_ad(V(:,i))'*convertVelocityToRegressor(V(:,i)));
        end

        % Y -> YB
        Y_B = Y*B'/(B*B');
        sum_A = sum_A + Y_B'*sigma_inv*Y_B;
        sum_A = (sum_A + sum_A')/2;
    end

    B_metric_invBt_half = sqrtm(B_metric_inv_Bt);
    B_metric_invBt_half = (B_metric_invBt_half+B_metric_invBt_half')/2;

    C = B_metric_invBt_half * sum_A * B_metric_invBt_half;
    C = (C+C')/2;

end