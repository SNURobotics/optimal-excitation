function [A_B, b] = getFullRegressorTorque(robot,p,trajectory)

    B = robot.B;
    num_base = size(B,1);
    
    cum_A = zeros(trajectory.num_sample*robot.dof,num_base);
    b              = zeros(trajectory.num_sample*robot.dof,1);

    sample_time = linspace(0,trajectory.horizon,trajectory.num_sample);

    % trajectory generation with parameter p
%    [q, qdot, qddot] = makeSpline(p, trajectory.order, trajectory.horizon, sample_time);
   [q, qdot, qddot] = makeFourier(p, trajectory.base_frequency, sample_time);

    for t=1:trajectory.num_sample
        [tau, V, Vdot] = solveInverseDynamics(robot.A,robot.M,q(:,t),qdot(:,t),qddot(:,t),robot.G);
        Y = getRegressorRecursive(robot.A,robot.M,q(:,t),V,Vdot);

        b(robot.dof*(t-1)+1:robot.dof*t) = tau;
        
        Y_B = Y*B'/(B*B');
        cum_A(robot.dof*(t-1)+1:robot.dof*t,:) = Y_B;
    end
    
    A_B = cum_A;
end