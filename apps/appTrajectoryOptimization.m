%% Trajectory Optimization for Parameter Identification

%% Initialization
n = 3;              % num of joints
m = 10;             % num of spline coefs per joints
basis_order = 4;    % cubic spline
total_time = 30;    % trajectory horizon


parameters = rand(m,n)
trajectory = makeSpline(parameters, basis_order, total_time);

figure()
hold on;
for t = 0:0.1:total_time
    plot(fnval(trajectory(1),t), fnval(trajectory(2),t), '.')
end