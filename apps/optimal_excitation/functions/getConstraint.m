%% Constraints and its Derivatives for Trajectory Optimization
% 2018 Bryan Dongik Lee

%% Inputs

%% Outputs

%% Implementation
function [c, ceq] = getConstraint(p, trajectory)
    ceq = [];
    
    % joint constraints
    q_max = ones(7,1)*pi';
    q_min = -q_max;
    qdot_max = ones(7,1)*pi';
    qdot_min = -qdot_max;
    qddot_max =  ones(7,1)*pi';
    qddot_min  = -qddot_max;
    
    % Fourier trajectory parameters
    num_sample       = trajectory.num_sample;
    horizon          = trajectory.horizon;
    base_frequency   = trajectory.base_frequency;
    sample_time      = linspace(0,horizon,num_sample);

    [q, qdot, qddot] = makeFourier(p, base_frequency, sample_time);
    c = [q-q_max, q_min-q, qdot-qdot_max, qdot_min-qdot, qddot-qddot_max, qddot_min-qddot];
end