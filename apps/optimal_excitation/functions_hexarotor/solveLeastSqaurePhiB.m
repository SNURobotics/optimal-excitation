function Phi_B = solveLeastSqaurePhiB(robot,p,trajectory, sigma)

    sigma_inv = pinv(sigma);
    Sigma_inv = zeros(trajectory.num_sample*robot.dof);
    for i = 1:trajectory.num_sample
        Sigma_inv(robot.dof*(i-1)+1:robot.dof*i,robot.dof*(i-1)+1:robot.dof*i) = sigma_inv;
    end

    B = robot.B;
    num_base = size(B,1);

    cummulativeY_B = zeros(trajectory.num_sample*robot.dof,num_base);
    b              = zeros(trajectory.num_sample*robot.dof,1);

    sample_time = linspace(0,trajectory.horizon,trajectory.num_sample);
    S_inv = robot.S_inv;

    % trajectory generation with parameter p
   [q, qdot, qddot] = makeSpline(p, trajectory.order, trajectory.horizon, sample_time);
%    [q, qdot, qddot] = makeFourier(p, trajectory.base_frequency, sample_time);

    for t=1:trajectory.num_sample
        [tau, V, Vdot] = solveInverseDynamics(robot.A,robot.M,q(:,t),qdot(:,t),qddot(:,t),robot.G);
        Y = S_inv*(convertVelocityToRegressor(Vdot(:,robot.dof)) - small_ad(V(:,robot.dof))'*convertVelocityToRegressor(V(:,robot.dof)));
        u = Y*robot.Phi;

        b(robot.dof*(t-1)+1:robot.dof*t) = u + sigma*randn(robot.dof,1);
        cummulativeY_B(robot.dof*(t-1)+1:robot.dof*t,:) = Y;
    end
    
%     cond_number = cond(cummulativeY_B'*Sigma_inv*cummulativeY_B*robot.B_metric_inv_Phi_Bt)
    
    Phi_B = pinv(cummulativeY_B'*Sigma_inv*cummulativeY_B)*cummulativeY_B'*Sigma_inv*b;
end