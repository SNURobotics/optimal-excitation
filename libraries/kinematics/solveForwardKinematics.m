%% Forward Kinematics Solver
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]     [Description]                       [Size]
%  q          joint angles                        n*1
%  A          i-th body screws from i-th frame    6*n
%  M          initial relative frames M_{i,i-1}   4*4*n

%% Outputs
% [Name]     [Description]                       [Size]
%  T          end-effector frame SE(3)            4*4

%% Examples

%% Implementation
function T = solveForwardKinematics(q,A,M)
    n = size(q,1); % number of joints
    T = eye(4,4);
    for i = 1:n
        T = T*inverse_SE3(M(:,:,i))*exp_se3(A(:,i)*q(i));
    end
end