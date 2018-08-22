%% SE(3) Visualizer
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]      [Description]                                 [Size]
%  T           frame SE(3) to draw                           4*4
%  frame       (optional) existing frame(handler) to update  struct created by this visualizer

%% Outputs
% [Name]      [Description]                                 [Size]
%  frame       frame handler for additional update           struct

%% Examples
% end_effector = plot_SE3(end_effector_T);
% plot_SE3(end_effector_T, end_effector);

%% Implementation
function frame = plot_SE3(T, varargin)
    if(size(T) ~= [4,4])
        error('se3 size error');
    end
    
    p=T(1:3,4);
    ax=p+T(1:3,1)*0.1;
    ay=p+T(1:3,2)*0.1;
    az=p+T(1:3,3)*0.1;
    
    if (nargin == 1)
        hold on;
        frame.x = plot3([p(1),ax(1)], [p(2),ax(2)], [p(3),ax(3)],'r');
        frame.y = plot3([p(1),ay(1)], [p(2),ay(2)], [p(3),ay(3)],'g');
        frame.z = plot3([p(1),az(1)], [p(2),az(2)], [p(3),az(3)],'b');
    elseif(nargin == 2)
        frame = varargin{1};
        set(frame.x, 'XData',[p(1),ax(1)], 'YData',[p(2),ax(2)], 'ZData',[p(3),ax(3)]);
        set(frame.y, 'XData',[p(1),ay(1)], 'YData',[p(2),ay(2)], 'ZData',[p(3),ay(3)]);
        set(frame.z, 'XData',[p(1),az(1)], 'YData',[p(2),az(2)], 'ZData',[p(3),az(3)]);
    end  
end