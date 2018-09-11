%% Sample a random pose within joint limits
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]      [Description]                                  [Size]
%  robot       robot object with joint limits                 struct

%% Outputs
% [Name]           [Description]                             [Size]
%  feasible_pose    random pose within joint limits           dof*1

%% Implementation
function feasible_pose = getFeasiblePose(robot)
    feasible_pose = zeros(robot.dof,1);
    for i = 1:robot.dof
        feasible_pose(i) = (rand() - 0.5)*(robot.q_max(i) - robot.q_min(i)) + (robot.q_max(i) + robot.q_min(i))/2;
    end
end