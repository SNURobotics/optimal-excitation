%% Compute the Center of Mass of a Tree Structure Robot
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]      [Description]                     [Size]
%  robot       robot object                      struct
%  q           joint angles                      n*1

%% Outputs
% [Name]           [Description]                                     [Size]
%  center_of_mass   center of mass of the robot from the base frame   3*1

%% Implementation
function center_of_mass = getCOM(robot, q)
    n = robot.dof;        % number of joints
    
    % kinematic variables
    T_global    = zeros(4,4,n); % T_i0 
    T_global(:,:,robot.root) = eye(4);

    m_total = 0;
    mp_total = zeros(3,1);
    
    stack = CStack();
    stack.push(robot.root);
    while(~stack.isempty)
        i = stack.pop;

        for child = 1:size(robot.tree{i}.children,1)
            stack.push(robot.tree{i}.children(child));
        end

        % T, Ad_T
        T_global(:,:,i) = T_global(:,:,robot.tree{i}.parent)*inverse_SE3(robot.M(:,:,i))*exp_se3(robot.A(:,i) * q(i));
        
        m = robot.Phi(1 + 10*(i-1));
        p = robot.Phi(2 + 10*(i-1):4 + 10*(i-1))/m;
        
        p_global = T_global(1:3,1:3,i)*p + T_global(1:3,4,i);
        
        m_total = m_total + m;
        mp_total = mp_total + m*p_global;
    end

    center_of_mass = mp_total/m_total;
end