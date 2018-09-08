%% Inertia Tensor Visualizer
% Created by TY Lee
% Modified by Bryan Lee

%% Inputs
% [Name]      [Description]                          [Size]
%  T_sb_i      body frames SE(3) from base frame      4*4*n
%  G_b_i       body inertia matrix G                  6*6*n        
%  trans       (optional) transparency                1*1
%  color       (optional) inertia color RGBs          n*3

%% Implementation
function plot_InertiaTensor(T_sb_i, G_b_i, varargin)
    n_parts = size(G_b_i,3);
    
    if nargin == 2
        trans = 0.2;
        color = rand(n_parts,3);
    elseif nargin == 4
        trans = varargin{1};
        color = varargin{2};
    end

    for i = 1:n_parts
    % [R,S,V] = svd(G_sbar_init{i}(1:3,1:3));
    % a = sqrt(5*(-S(1,1)+S(2,2)+S(3,3))/(2*G_sbar_init{i}(4,4)));
    % b = sqrt(5*(S(1,1)-S(2,2)+S(3,3))/(2*G_sbar_init{i}(4,4)));
    % c = sqrt(5*(S(1,1)+S(2,2)-S(3,3))/(2*G_sbar_init{i}(4,4)));
    T_bc = [eye(3,3), skew(G_b_i(1:3,4:6,i))./G_b_i(4,4,i);0,0,0,1];
    G_c = large_Ad(T_bc)'*G_b_i(:,:,i)*large_Ad(T_bc);
    % [R,S,V] = svd(G_b(1:3,1:3));
    [R,S,V] = eig(G_c(1:3,1:3));

    if det(R)<0
        R_temp = R; S_temp = S;
        R(:,1) = R_temp(:,2);
        R(:,2) = R_temp(:,1);
        S(1,1) = S_temp(2,2);
        S(2,2) = S_temp(1,1);
    end
    T = [R, zeros(3,1);0,0,0,1];
    T(1:3,4) = skew(G_b_i(1:3,4:6,i))./G_b_i(4,4,i);
    T = T_sb_i(:,:,i)*T;
    a = sqrt(5*(-S(1,1)+S(2,2)+S(3,3))/(2*(G_b_i(4,4,i))));
    if ~isreal(a)
        a = 0;
    end
    b = sqrt(5*(S(1,1)-S(2,2)+S(3,3))/(2*(G_b_i(4,4,i))));
    if ~isreal(b)
        b = 0;
    end
    c = sqrt(5*(S(1,1)+S(2,2)-S(3,3))/(2*(G_b_i(4,4,i))));
    if ~isreal(c)
        c = 0;
    end
    if a*b*c ~= 0
    % T(1:3,1:3) = T(1:3,1:3)*[a,0,0;0,b,0;0,0,c];
    % draw_SE3(T);
    % transl = T(1:3,1:3)'*T(1:3,4);
    % density = G_b_i(4,4,i)/((a+0.001)*(b+0.001)*(c+0.001)*2e3);

    [xc,yc,zc] = ellipsoid(0,0,0,a,b,c);

    % rotate data with orientation matrix U and center M
    RR = T(1:3,1:3); M = T(1:3,4);
    a = kron(RR(:,1),xc); b = kron(RR(:,2),yc); c = kron(RR(:,3),zc);
    data = a+b+c; n = size(data,2);
    x = data(1:n,:)+M(1); y = data(n+1:2*n,:)+M(2); z = data(2*n+1:end,:)+M(3);
    % scatter3(M(1),M(2),M(3), 'MarkerEdgeColor',color(i,:),'MarkerFaceColor',color(i,:));
    % now plot the rotated ellipse
    % color_z = ones(size(z))*color;


    sc = surf(x,y,z,'EdgeColor','none','FaceColor',color(i,:),'FaceAlpha',trans);
    else
        M = T(1:3,4);
        text(M(1),M(2),M(3),['     Link' num2str(i) 'infeasible'],'HorizontalAlignment','left','FontSize',8);
    end
    hold on;
    end
    axis equal; hold off;
end