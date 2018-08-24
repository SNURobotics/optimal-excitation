close all
clear
clc

robot = makeHexarotor();

q = zeros(robot.dof,1);
T = solveForwardKinematics(q, robot.A, robot.M);

%% Figure Setting
figure('Name','Phantom','NumberTitle','off','units','pixels','pos',[100 100 1000 1000]);
hold on;
axis equal;
axis([-2 2 -2 2 -2 2]);
xlabel('x'); ylabel('y'); zlabel('z');
view([-45 35]);

%% STL Load
fv_zero = stlread('Phantom.STL');  % base link
fv_zero.vertices = (T(1:3,1:3)*fv_zero.vertices' + T(1:3,4)*ones(1,size(fv_zero.vertices,1)))';

%% STL 
fv = fv_zero;
render = patch(fv,'FaceColor',       [0.7 0.7 0.7], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15);
     
camlight('headlight');
material('dull');

uav = plot_SE3(T);

plot_InertiaTensor(eye(4), robot.G);

getframe;

%% Animation
m = 10;
p = 2*rand(m,robot.dof);

num_sample = 200;
sample_time      = linspace(0,20,num_sample);

q = makeSpline(p, 4, 20, sample_time);

while(1)
    for i = 1:num_sample
        T = solveForwardKinematics(q(:,i), robot.A, robot.M);

        fv.vertices = (T(1:3,1:3)*fv_zero.vertices' + T(1:3,4)*ones(1,size(fv_zero.vertices,1)))';
        set(render, 'Vertices', fv.vertices);

        plot_SE3(T, uav);
        getframe;
    end
end