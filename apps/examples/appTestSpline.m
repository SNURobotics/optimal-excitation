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
trajectory  = makeSpline(parameters, basis_order, total_time);
p2ptrajectory = makeSplineP2P(qi,qf,parameters, basis_order, total_time);

trajectory_derivative = getSplineDerivative(parameters, basis_order, total_time);

%%
tr = zeros(n,10);
tr_from_derivative = zeros(n,10);
for t= 1:10
    for i = 1:n
        tr(i,t) = fnval(trajectory(i), t);    

        for j = 1:m
            tr_from_derivative(i,t) = tr_from_derivative(i,t) + fnval(trajectory_derivative(i,j,i),t)*parameters(j,i);
        end
    end
end 

tr
tr_from_derivative

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

%% 5th Order Spline
qi = rand(4,1);
qidot = rand(4,1);
qiddot = rand(4,1);
qf = rand(4,1);
qfdot = rand(4,1);
qfddot = rand(4,1);

horizon = 10;
sample_time = linspace(0,horizon,100);
[q, qdot, qddot] = make5thSplineP2P(qi, qidot, qiddot, qf, qfdot, qfddot, horizon, sample_time);

figure();
plot(sample_time, q); hold on;
plot(0, qi, '*');
plot(horizon, qf, '*');

plot(sample_time, qdot); hold on;
plot(0, qidot, '*');
plot(horizon, qfdot, '*');

plot(sample_time, qddot); hold on;
plot(0, qiddot, '*');
plot(horizon, qfddot, '*');
