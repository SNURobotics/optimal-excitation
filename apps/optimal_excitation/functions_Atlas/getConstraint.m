%% Constraints and its Derivatives for Trajectory Optimization
% 2018 Bryan Dongik Lee

%% Inputs

%% Outputs

%% Implementation
function [c, ceq, gradc, gradceq] = getConstraint(p, trajectory, robot)
    ceq = [];
    
    % joint constraints
    q_max = robot.q_max;
    q_min = robot.q_min;
    qdot_max = robot.qdot_max;
    qdot_min = robot.qdot_min;
    
    % Fourier trajectory parameters
    num_sample       = trajectory.num_sample;
    horizon          = trajectory.horizon;
    base_frequency   = trajectory.base_frequency;
    sample_time      = linspace(0,horizon,num_sample);

    [q, qdot] = makeFourier(p, base_frequency, sample_time, (robot.q_max + robot.q_min)/2);
    c = [q-q_max, q_min-q, qdot-qdot_max, qdot_min-qdot];
    
    if nargout > 2
        gradceq = [];
        
        m = size(p,1);
        n = robot.dof;
        
        [dq, dqdot] = getFourierDerivative(p, base_frequency, sample_time);
        
        gradc = zeros(robot.dof, 4*num_sample, m, n);
        for p = 1:m
            for l = 1:n
                gradc(:,:,p,l) = [dq(:,p,l,:) -dq(:,p,l,:) dqdot(:,p,l,:) -dqdot(:,p,l,:)];
            end
        end
    end
end