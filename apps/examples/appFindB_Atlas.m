close all
clear 
clc

robot = makeAtlasV5();

n = robot.dof;
A = robot.A;
M = robot.M;

num_sample = 3000;

T           = zeros(4,4,n); % T_{i,i-1}
Ad_T        = zeros(6,6,n); % Ad_T_{i,i-1}
T_global    = zeros(4,4,n); % T_0i 
T_global(:,:,robot.root) = eye(4);
Ad_T_global = zeros(6,6,n); % Ad_T_i0
Ad_T_global(:,:,robot.root) = eye(6);

V_0 = zeros(6,1);
Vdot_0 = zeros(6,1);
Vdot_0(6) = 9.8;

V = zeros(6,n);
Vdot = zeros(6,n);

Y = zeros(6*num_sample,10*n);

for iter = 1:num_sample
    
    q = (rand(n,1) - 0.5)*pi*2;
    qdot = (rand(n,1) - 0.5)*pi*100;
    qddot = (rand(n,1) - 0.5)*pi*100;

    stack = CStack();
    stack.push(robot.root);
    while(~stack.isempty)
        i = stack.pop;
        
        for child = 1:size(robot.tree{i}.children,1)
            stack.push(robot.tree{i}.children(child));
        end

        T(:,:,i) = exp_se3(-A(:,i)*q(i))*M(:,:,i);
        Ad_T(:,:,i) = large_Ad(T(:,:,i));
        T_global(:,:,i) = T_global(:,:,robot.tree{i}.parent)*inverse_SE3(robot.M(:,:,i))*exp_se3(robot.A(:,i) * q(i));
        Ad_T_global(:,:,i) = large_Ad(inverse_SE3(T_global(:,:,i)));

        if i == robot.root
            V(:,i)    = V_0;
            Vdot(:,i) = Vdot_0;
        else
            V(:,i)    = Ad_T(:,:,i)*V(:,robot.tree{i}.parent)    + A(:,i)*qdot(i);
            Vdot(:,i) = Ad_T(:,:,i)*Vdot(:,robot.tree{i}.parent) + small_ad(V(:,i))*A(:,i)*qdot(i) + A(:,i)*qddot(i);
        end
    end
    
    for i = 1:n
        Y(6*(iter-1)+1:6*iter, 10*(i-1)+1:10*i) = Ad_T_global(:,:,i)'...
                    *(convertVelocityToRegressor(Vdot(:,i)) - small_ad(V(:,i))'*convertVelocityToRegressor(V(:,i)));
    end
end

[U,S,V] = svd(Y);
B = V(:,1:204)'