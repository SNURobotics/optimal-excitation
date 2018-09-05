close all
clear
clc

robot = makeAtlasV5();

dof = robot.dof;
T = zeros(4,4,dof);
T(:,:,1) = eye(4);
q = zeros(dof,1);


for i = 1:dof
    if strcmp(robot.name{i}, 'r_lfarm')
        robot.stl{i} = stlread('r_hand.stl');
    elseif strcmp(robot.name{i}, 'r_ufarm')
        robot.stl{i} = stlread('r_farm.stl');
    elseif strcmp(robot.name{i}, 'r_hand')
        robot.stl{i} = stlread('r_hand.stl');

    elseif strcmp(robot.name{i}, 'l_lfarm')
        robot.stl{i} = stlread('r_hand.stl');
    elseif strcmp(robot.name{i}, 'l_ufarm')
        robot.stl{i} = stlread('r_farm.stl');
    elseif strcmp(robot.name{i}, 'l_hand')
        robot.stl{i} = stlread('r_hand.stl');
        
    elseif strcmp(robot.name{i}, 'l_clav')
        robot.stl{i} = stlread('r_clav.stl');
    elseif strcmp(robot.name{i}, 'l_scap')
        robot.stl{i} = stlread('r_scap.stl');
    elseif strcmp(robot.name{i}, 'l_uarm')
        robot.stl{i} = stlread('r_uarm.stl');
    elseif strcmp(robot.name{i}, 'l_larm')
        robot.stl{i} = stlread('r_larm.stl');        
    else
        robot.stl{i} = stlread([robot.name{i} '.stl']);
    end
end

figure('Name','Atlas v5','NumberTitle','off','units','pixels','pos',[-1000 200 900 900]);
hold on;
axis equal;
axis([-2 2 -2 2 -1 3]);
xlabel('x'); ylabel('y'); zlabel('z');
view(135, 30);
         
% Add a camera light, and tone down the specular highlighting
camlight('headlight');
material('dull');

stack = CStack();
stack.push(1);

while(~stack.isempty)
    link = stack.pop;
    
    for i = 1:size(robot.tree{link}.children,1)
        stack.push(robot.tree{link}.children(i));
    end

    try
        T(:,:,link) = T(:,:,robot.tree{link}.parent)*inverse_SE3(robot.M(:,:,link))*exp_se3(robot.A(:,link) * q(link));
    catch 
        if robot.tree{link}.parent == 0
            T(:,:,link) = inverse_SE3(robot.M(:,:,link))*exp_se3(robot.A(:,link) * q(link));
        end
    end
    
    plot_SE3(T(:,:,link));
    
    robot.stl{link}.vertices = (T(1:3,1:3,link)*robot.stl{link}.vertices' ...
                                + T(1:3,4,link)*ones(1,size(robot.stl{link}.vertices,1)))';

    render_part{link} = patch(robot.stl{link},'FaceColor',  [0.7 0.7 0.7], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15);
     
     plot_inertiatensor(T(:,:,link), robot.G(:,:,link));
    
end