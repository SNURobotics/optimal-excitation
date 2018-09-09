%% Visualizer for Atlas V5
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]       [Description]                    [Size]
%  trajectory   robot trajectory to visualize    7*num_time

%% Implementation
function appVisualizeAtlas(varargin)    

    robot = makeAtlasV5();

    dof = robot.dof;
    T = zeros(4,4,dof);
    T(:,:,robot.root) = eye(4);

    trajectory = zeros(dof,1);
    if nargin > 1
        if size(varargin{1},1) == dof
            trajectory = varargin{1};
        else
            error('number of the joints and the trajectory does not match');
        end
    end   
    
    num_time = size(trajectory, 2);

   
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
    axis([-1.5 1.5 -1.5 1.5 -0.5 2.5]);
    xlabel('x'); ylabel('y'); zlabel('z');
    view(110, 30);

    % Add a camera light, and tone down the specular highlighting
    camlight('headlight');
    material('dull');

    stack = CStack();
    stack.push(robot.root);
    while(~stack.isempty)
        link = stack.pop;

        for i = 1:size(robot.tree{link}.children,1)
            stack.push(robot.tree{link}.children(i));
        end

        T(:,:,link) = T(:,:,robot.tree{link}.parent)*inverse_SE3(robot.M(:,:,link))*exp_se3(robot.A(:,link) * trajectory(link,1));

        plot_SE3(T(:,:,link));
        robot.stl_zero{link} = robot.stl{link};
        robot.stl{link}.vertices = (T(1:3,1:3,link)*robot.stl{link}.vertices' ...
                                    + T(1:3,4,link)*ones(1,size(robot.stl{link}.vertices,1)))';

        render_part{link} = patch(robot.stl{link},'FaceColor',  [0.7 0.7 0.7], ...
             'EdgeColor',       'none',        ...
             'FaceLighting',    'gouraud',     ...
             'AmbientStrength', 0.15);

%          plot_inertiatensor(T(:,:,link), robot.G(:,:,link));   
    end


    while(1)
        for time = 1:num_time
            stack.push(robot.root);
            while(~stack.isempty)
                link = stack.pop;

                for i = 1:size(robot.tree{link}.children,1)
                    stack.push(robot.tree{link}.children(i));
                end

                T(:,:,link) = T(:,:,robot.tree{link}.parent)*inverse_SE3(robot.M(:,:,link))*exp_se3(robot.A(:,link) * trajectory(link,time));
                robot.stl{link}.vertices = (T(1:3,1:3,link)*robot.stl_zero{link}.vertices' ...
                                            + T(1:3,4,link)*ones(1,size(robot.stl_zero{link}.vertices,1)))';
                set(render_part{link}, 'Vertices', robot.stl{link}.vertices);
            end
            getframe;
        end
    end
end