%% Sample random points within joint limits

function feasible_initial = getInit(robot)
    feasible_initial = zeros(robot.dof,1);
    for i = 1:robot.dof
        feasible_initial(i) = (rand() - 0.5)*(robot.q_max(i) - robot.q_min(i)) + (robot.q_max(i) + robot.q_min(i))/2;
    end
end