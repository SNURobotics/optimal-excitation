%% Forward Kinematics Solver
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]       [Description]                      [Size]

%% Outputs
% [Name]  [Description]                      [Size]

%% Examples

%% Implementation
function tool_T = solveForwardKinematics(q,A,M)
    
    n = size(q,1); % number of joints
    tool_T = eye(4,4);
    for i = 1:n
        tool_T = tool_T*inverse_SE3(M(:,:,i))*exp_se3(A(:,i)*q(i));
    end
    
end