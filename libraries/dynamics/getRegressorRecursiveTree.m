%% Recursive Regressor Calculator for Tree Structure
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]       [Description]                      [Size]
%  robot        robot object w/ tree               struct
%  Ad_T         Ad_T{i, i-1}                       6*6*n
%  V            spatial velocities                 6*n
%  Vdot         spatial accelerations              6*n

%% Outputs
% [Name]  [Description]                                     [Size]
%  Y       regressor matrix (tau = Y*PHI)                    n*10n

%% Implementation
function Y = getRegressorRecursiveTree(robot, V, Vdot, Ad_T)
    n= robot.dof;
    W = zeros(6*n,10*n);

    diagA  = zeros(6*n,n);
    for i=1:n
        diagA(6*i-5:6*i,i) = robot.A(:,i);
    end 
    
    %%
    leaf = [];
    queue = CQueue();
    stack = CStack();
    queue.push(robot.root);

    while ~queue.isempty
        i = queue.pop();
        stack.push(i);

        if size(robot.tree{i}.children,1) > 0
            for child = 1:size(robot.tree{i}.children,1)
                    queue.push(robot.tree{i}.children(child));
            end
        else
            leaf = [leaf; i];
        end
    end

    while ~stack.isempty
        i = stack.pop();
            
        W(6*i-5:6*i,10*i-9:10*i) = convertVelocityToRegressor(Vdot(:,i)) - small_ad(V(:,i))'*convertVelocityToRegressor(V(:,i));
        
        for j = 1:size(robot.tree{i}.children,1)
            stack_local = CStack();
            stack_local.push(robot.tree{i}.children(j));
            
            while(~stack_local.isempty)
                k = stack_local.pop;

                for child = 1:size(robot.tree{k}.children,1)
                    stack_local.push(robot.tree{k}.children(child));
                end
                
                W(6*i-5:6*i,10*k-9:10*k) = Ad_T(:,:,robot.tree{i}.children(j))'*W(6*robot.tree{i}.children(j)-5:6*robot.tree{i}.children(j),10*k-9:10*k);                
            end       
        end
    end

    Y = diagA'*W;
end