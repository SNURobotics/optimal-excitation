%% Compute the Center of Mass of a Tree Structure Robot
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]      [Description]                     [Size]
%  robot       robot object (A, B, M)            struct
%  q           joint angles                      n*1

%% Outputs
% [Name]    [Description]           [Size]
%  zmp       zero momentum point     2*1

%% Implementation
function zmp = getZMP(robot, q, qdot, qddot)
    n = robot.dof;        % number of joints
    
    % kinematic variables
    T           = zeros(4,4,n); % T_{i,i-1}
    Ad_T        = zeros(6,6,n); % Ad_T_{i,i-1}
    T_global    = zeros(4,4,n); % T_i0 
    T_global(:,:,robot.root) = eye(4);
    Ad_T_global = zeros(6,6,n); % Ad_T_i0
    Ad_T_global(:,:,robot.root) = eye(6);
    
    V_0     = zeros(6,1);       % base velocity
    Vdot_0  = zeros(6,1);       % base acceleration
    Vdot_0(6) = 9.8;
    
    V = zeros(6,n);
    Vdot = zeros(6,n);
    
    F = zeros(6,1);
    
    stack = CStack();
    stack.push(robot.root);
    while(~stack.isempty)
        i = stack.pop;

        for child = 1:size(robot.tree{i}.children,1)
            stack.push(robot.tree{i}.children(child));
        end

        % T, Ad_T
        T(:,:,i) = exp_se3(-robot.A(:,i)*q(i))*robot.M(:,:,i);
        Ad_T(:,:,i) = large_Ad(T(:,:,i));
        T_global(:,:,i) = exp_se3(-robot.A(:,i) * q(i))* robot.M(:,:,i) *T_global(:,:,robot.tree{i}.parent);
        Ad_T_global(:,:,i) = large_Ad(T_global(:,:,i));
            
        % V, Vdot
        if i == robot.root
            V(:,i)    = V_0;
            Vdot(:,i) = Vdot_0;
        else
            V(:,i)    = Ad_T(:,:,i)*V(:,robot.tree{i}.parent)    + robot.A(:,i)*qdot(i);
            Vdot(:,i) = Ad_T(:,:,i)*Vdot(:,robot.tree{i}.parent) + small_ad(V(:,i))*robot.A(:,i)*qdot(i) + robot.A(:,i)*qddot(i);
        end
        
        
        F = F + Ad_T_global(:,:,i)'*(robot.G(:,:,i)*Vdot(:,i) - small_ad(V(:,i))'*robot.G(:,:,i)*V(:,i));
    end
    
    zmp = [-F(2)/F(6); F(1)/F(6)];
end