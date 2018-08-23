close all
clear
clc

%% Robot Init
robot = makeKukaR820();

n = robot.dof;
q = zeros(n,1);
T = zeros(4,4,n);
for i = 1:n
    T(:,:,i) = solveForwardKinematics(q(1:i), robot.A(:,1:i), robot.M(:,:,1:i));
end

%% STL Load
fv_base = stlread(['link_0.STL']);  % base link
fv_base.vertices = (T(1:3,1:3,1)*fv_base.vertices' + T(1:3,4,1)*ones(1,size(fv_base.vertices,1)))';

for i = 1:n
    fv_zero{i} = stlread(['link_' num2str(i) '.STL']);
end

fv = fv_zero;
for i = 1:n
    fv{i}.vertices = (T(1:3,1:3,i)*fv{i}.vertices' + T(1:3,4,i)*ones(1,size(fv{i}.vertices,1)))';
end

%% Render
% The model is rendered with a PATCH graphics object. We also add some dynamic
% lighting, and adjust the material properties to change the specular
% highlighting.
figure('Name','KUKA LWR iiwa R820','NumberTitle','off','units','pixels','pos',[100 100 1000 1000]);
hold on;
axis equal;
axis([-1 1 -1 1 -0.5 1.2]);
xlabel('x'); ylabel('y'); zlabel('z');

% draw base link
patch(fv_base,'FaceColor',       [0.7 0.7 0.7], ...
             'EdgeColor',       'none',        ...
             'FaceLighting',    'gouraud',     ...
             'AmbientStrength', 0.15);

% draw 7 links
color = ones(n,3)*0.8;
color(2,:) = [246 120 40]/255;
color(6,:) = [246 120 40]/255;

for i = 1:7
    render_part{i} = patch(fv{i},'FaceColor',  color(i,:), ...
             'EdgeColor',       'none',        ...
             'FaceLighting',    'gouraud',     ...
             'AmbientStrength', 0.15);
end

% draw end-effector
end_effector_M = eye(4);
end_effector_M(3,4) = 0.125;
end_effector_T = T(:,:,7) * end_effector_M;
end_effector = plot_SE3(end_effector_T);

% Add a camera light, and tone down the specular highlighting
camlight('headlight');
material('dull');

% Fix the axes scaling, and set a nice view angle
view([-135 35]);
getframe;

%% Animation
joint  = 1; 

while(1)
    q(joint) = q(joint) + 0.2;
    
    for i = 1:n
        T(:,:,i) = solveForwardKinematics(q(1:i), robot.A(:,1:i), robot.M(:,:,1:i));
        fv{i}.vertices = (T(1:3,1:3,i)*fv_zero{i}.vertices' + T(1:3,4,i)*ones(1,size(fv_zero{i}.vertices,1)))';
    end
    
    for i = 1:n
        set(render_part{i}, 'Vertices', fv{i}.vertices);
    end
    
    end_effector_T = T(:,:,7) * end_effector_M;
    plot_SE3(end_effector_T, end_effector);
    
    getframe;
end