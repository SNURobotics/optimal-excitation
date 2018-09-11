function Phi_B = solveLeastSqaurePhiB(robot,p,trajectory,sigma,max_variance)

    sigma_inv = pinv(sigma);
    Sigma_inv = zeros(trajectory.num_sample*robot.dof);
    for i = 1:trajectory.num_sample
        Sigma_inv(robot.dof*(i-1)+1:robot.dof*i,robot.dof*(i-1)+1:robot.dof*i) = sigma_inv;
    end

    B = robot.B;
    num_base = size(B,1);
    
    B_metric_inv_Bt = robot.B_metric_inv_Phi_Bt_0;
    
    sum_A_sigma_At = zeros(num_base, num_base);
    cum_A = zeros(trajectory.num_sample*robot.dof,num_base);
    b              = zeros(trajectory.num_sample*robot.dof,1);

    sample_time = linspace(0,trajectory.horizon,trajectory.num_sample);

    % trajectory generation with parameter p
%    [q, qdot, qddot] = makeSpline(p, trajectory.order, trajectory.horizon, sample_time);
   [q, qdot, qddot] = makeFourier(p, trajectory.base_frequency, sample_time);

    for t=1:trajectory.num_sample
        [tau, V, Vdot] = solveInverseDynamics(robot.A,robot.M,q(:,t),qdot(:,t),qddot(:,t),robot.G);
        Y = getRegressorRecursive(robot.A,robot.M,q(:,t),V,Vdot);

        b(robot.dof*(t-1)+1:robot.dof*t) = tau + sqrtm(sigma)*randn(robot.dof,1);
        
        Y_B = Y*B'/(B*B');
        cum_A(robot.dof*(t-1)+1:robot.dof*t,:) = Y_B;
        sum_A_sigma_At = sum_A_sigma_At + Y_B' * sigma_inv * Y_B;
    end
    
    B_metric_invBt_half = sqrtm(B_metric_inv_Bt);
    B_metric_invBt_half = (B_metric_invBt_half+B_metric_invBt_half')/2;

    C = B_metric_invBt_half* sum_A_sigma_At * B_metric_invBt_half;
    C = (C+C')/2;    
    
    % eigen decomposition
    [Q, D] = eig(C);
    D = inv(D);

    k = size(D,1);
    eigenvalues = zeros(k,1);
    for i=1:k
        eigenvalues(i) = D(i,i);
    end

    V = B_metric_invBt_half * Q;

    num_identify = size(find(eigenvalues <= max_variance),1);

    V_neg = zeros(size(V,1), num_identify);
    V_pos = zeros(size(V,1), size(V,1) - num_identify);
    lambda_neg = zeros(num_identify); 
    lambda_pos = zeros(size(V,1) - num_identify);

    neg = 0;
    pos = 0;
    for i = 1:size(V,1)
        if D(i,i) < max_variance
            neg = neg + 1;
            V_neg(:,neg) = V(:,i);
            lambda_neg(neg,neg) = D(i,i);
        else
            pos = pos + 1;
            V_pos(:,pos) = V(:,i);
            lambda_pos(pos,pos) = D(i,i);
        end
    end
    
    Phi_B = V_neg*lambda_neg*V_neg'*cum_A'*Sigma_inv*b + V_pos*V_pos'*pinv(B_metric_inv_Bt)*B*robot.Phi_0;    
end