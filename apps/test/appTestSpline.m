close all
clear
clc

n           = 2;    % num of joints
m           = 10;   % num of spline coefs per joints
basis_order = 4;    % cubic spline
total_time  = 30;   % trajectory horizon

qi = rand(n,1);
qf = rand(n,1);
parameters = rand(m,n);
trajectory  = makeSpline(parameters, basis_order, total_time, [1,2,3,4,5,6,7,8,9,10]);
p2ptrajectory = makeSplineP2P(qi,qf,parameters, basis_order, total_time);

trajectory_derivative = getSplineDerivative(parameters, basis_order, total_time, [1,2,3,4,5,6,7,8,9,10]);
tr = zeros(n,10);
for t= 1:10
    for i = 1:n
        for j = 1:m
            tr(i,t) = tr(i,t) + trajectory_derivative(i,j,t)*parameters(j,i);
        end
    end
end 

trajectory
tr

%%
% trajectory plot
figure(1)
title('Parmetric Spline')
hold on;
for t = 0:0.1:total_time
    plot(fnval(trajectory(1),t), fnval(trajectory(2),t), '.')
end

figure(2)
title('Point-to-Point Parametric Spline')
hold on;
plot(qi(1),qi(2),'*');
plot(qf(1),qf(2),'*');
for t = 0:0.1:total_time
    plot(fnval(p2ptrajectory(1),t), fnval(p2ptrajectory(2),t), '.')
end