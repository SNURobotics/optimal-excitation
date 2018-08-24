%% Visualizer for Hexarotor 
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]       [Description]                    [Size]
%  trajectory   robot trajectory to visualize    6*num_time

%% Implementation
function appVisualizeHexarotor(varargin)
    %% Robot Init
    robot = makeHexarotor();

    n = robot.dof;
    q = zeros(robot.dof,1);
    T = solveForwardKinematics(q, robot.A, robot.M);
    T_offset = eye(4); T_offset(3,4) = 0.2;
    T_coor = T*T_offset;
    
    trajectory = zeros(6,1);
    if nargin > 0
        if size(varargin{1},1) == n
            trajectory = varargin{1};
        else
            error('number of the joints and the trajectory does not match');
        end
    end
    
    num_time = size(trajectory, 2);
    
    %% Figure Setting
    figure('Name','Phantom','NumberTitle','off','units','pixels','pos',[100 100 1000 1000]);
    hold on;
    axis equal;
    axis([-1 1 -1 1 -1 1]*5);
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

    uav = plot_SE3(T_coor);

    getframe;

    %% Animation
    while(1)
        for time = 1:num_time
            T = solveForwardKinematics(trajectory(:,time), robot.A, robot.M);

            fv.vertices = (T(1:3,1:3)*fv_zero.vertices' + T(1:3,4)*ones(1,size(fv_zero.vertices,1)))';
            set(render, 'Vertices', fv.vertices);

            T_coor = T*T_offset;
            plot_SE3(T_coor, uav);
            plot3(T(1,4), T(2,4), T(3,4),'.b');
            getframe;
        end
    end
end