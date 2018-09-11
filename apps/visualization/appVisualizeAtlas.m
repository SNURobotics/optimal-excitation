%% Visualizer for Atlas V5
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]       [Description]                               [Size]
%  trajectory   (optional) robot trajectory to visualize    dof*num_time
%  speed        (optional) play speed                       1*1

%% Implementation
function appVisualizeAtlas(varargin)    
    %% options
    export_video = false;
    
    go_crazy = false;
    crazy = 1;
    
    
    %%    
    robot = makeAtlasV5();

    dof = robot.dof;
    T = zeros(4,4,dof);
    T(:,:,robot.root) = eye(4);

    % default pose and speed
    trajectory = zeros(dof,1);
    speed = 1;
    
    if nargin > 0
        if size(varargin{1},1) == dof
            trajectory = varargin{1};
        else
            error('number of the joints and the trajectory does not match');
        end
    
        if nargin > 1
            speed = varargin{2};
        end
    end   
    
    num_time = size(trajectory, 2);

    % load stl
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

    % window setting
    window = figure('Name','Atlas v5','NumberTitle','off','units','pixels','pos',[-1000 200 900 900]);
    hold on;
    axis equal;
    axis([-1.5 1.5 -1.5 1.5 -0.5 2.5]);
    axis off;
    xlabel('x'); ylabel('y'); zlabel('z');
    view(95, 10); % view angle
    bg_color = [0.6 0.8 1]; % figure background color
    set(gcf, 'color', bg_color);
    
    % Add a camera light, and tone down the specular highlighting
    camlight('headlight');
    material('dull');

    % color parts
    link_color = ones(robot.dof,3) * 0.9;
    link_color(1,:) = [0.5 0.5 0.5];
    link_color(8,:) = [1 0 0];
    link_color(2,:) = [1 0 0];
    link_color(12,:) = [1 0 0];
    link_color(6,:) = [1 0 0];
    
    % video export
    if export_video
        writerObj = VideoWriter('atlas','MPEG-4'); % 
        writerObj.FrameRate = 15;
    
        % open the video writer
        open(writerObj);
    end

    % initial rendering
    stack = CStack();
    stack.push(robot.root);
    while(~stack.isempty)
        link = stack.pop;

        for i = 1:size(robot.tree{link}.children,1)
            stack.push(robot.tree{link}.children(i));
        end

        T(:,:,link) = T(:,:,robot.tree{link}.parent)*inverse_SE3(robot.M(:,:,link))*exp_se3(robot.A(:,link) * trajectory(link,1));

%         plot_SE3(T(:,:,link));
        robot.stl_zero{link} = robot.stl{link};
        robot.stl{link}.vertices = (T(1:3,1:3,link)*robot.stl{link}.vertices' ...
                                    + T(1:3,4,link)*ones(1,size(robot.stl{link}.vertices,1)))';

        render_part{link} = patch(robot.stl{link},'FaceColor',  link_color(link,:), ...
             'EdgeColor',       'none',        ...
             'FaceLighting',    'gouraud',     ...
             'AmbientStrength', 0.15,          ...
             'FaceAlpha', 0.7);

%          plot_inertiatensor(T(:,:,link), robot.G(:,:,link));
    end
    
    % animation
    for i = 1:1
        for time = 1:speed:num_time
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
            frame = getframe(gcf);
            
            if export_video
                writeVideo(writerObj, frame);
            end
            
            if go_crazy && crazy > 1
                crazy = 1;
                bg_color = bg_color + randn(1,3);
                bg_color = bg_color - floor(bg_color);
                set(gcf, 'color', bg_color) % figure background color
            end
            crazy =  crazy + 1;
        end
    end
    
    if export_video
        % close the writer object
        close(writerObj);
    end
end