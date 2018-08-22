close all
clear
clc

robot = makeKukaR820;

figure('units','pixels','pos',[-1000 200 900 900]); hold on;
axis([-1 1 -1 1 -1 1]);
xlabel('x'); ylabel('y'); zlabel('z');
view(-45, 30);

base = plot_SE3(eye(4));

q = zeros(7,1);
end_effector_T = solveForwardKinematics(q, robot.A, robot.M);
end_effector = plot_SE3(end_effector_T);

getframe;

%%
joint = 1;
while(1)
    q(joint) = q(joint) + 0.1;
    end_effector_T = solveForwardKinematics(q, robot.A, robot.M);

    plot_SE3(end_effector_T, end_effector);
    getframe;
    pause(0.01);
end