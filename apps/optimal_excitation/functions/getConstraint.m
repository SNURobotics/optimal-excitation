%% Constraints and its Derivatives for Trajectory Optimization
% 2018 Bryan Dongik Lee

%% Inputs

%% Outputs

%% Implementation
function [c, ceq] = getConstraint(p, trajectory, robot)
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

    [q, qdot] = makeFourier(p, base_frequency, sample_time);
    c = [q-q_max, q_min-q, qdot-qdot_max, qdot_min-qdot];
end