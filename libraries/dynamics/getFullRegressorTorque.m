function [A_B, b] = getFullRegressorTorque(robot,p,trajectory)

    B = robot.B;
    num_base = size(B,1);
    
    cum_A = zeros(trajectory.num_sample*robot.dof,num_base);
    b              = zeros(trajectory.num_sample*robot.dof,1);

    sample_time = linspace(0,trajectory.horizon,trajectory.num_sample);

    % trajectory generation with parameter p
%    [q, qdot, qddot] = makeSpline(p, trajectory.order, trajectory.horizon, sample_time);
   [q, qdot, qddot] = makeFourier(p, trajectory.base_frequency, sample_time);
   
   %%%
   MassMat = zeros(robot.dof, robot.dof);
   %%%
    for t=1:trajectory.num_sample
        [tau, V, Vdot] = solveInverseDynamics(robot.A,robot.M,q(:,t),qdot(:,t),qddot(:,t),robot.G, [0;0;0;0;0;9.8]);
        Y = getRegressorRecursive(robot.A,robot.M,q(:,t),V,Vdot);
        
        b(robot.dof*(t-1)+1:robot.dof*t) = tau;
        
        Y_B = Y*B'/(B*B');
        
        %%%
        for k = 1 : robot.dof
            q_k = zeros(robot.dof,1); q_k(k) = 1;
           [MassMat_col,~,~] = solveInverseDynamics(robot.A,robot.M,q(:,t),zeros(robot.dof,1),q_k, robot.G, [0;0;0;0;0;9.8]);
           MassMat(:,k) = MassMat_col;
        end
        MassMat = MassMat^(0.5);
        Y_B = inv(MassMat)*Y_B;
        b(robot.dof*(t-1)+1:robot.dof*t) = inv(MassMat)*b(robot.dof*(t-1)+1:robot.dof*t);
%         %%%
        
        cum_A(robot.dof*(t-1)+1:robot.dof*t,:) = Y_B;
    end
    
    A_B = cum_A;
end