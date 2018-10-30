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
function M = getInertiaMatrixMTree(robot, Ad_T)
    n= robot.dof;
    
    largeA  = zeros(6*n,n);
    for i=1:n
        largeA(6*i-5:6*i,i) = robot.A(:,i);
    end 
    
    largeG = zeros(6*n,6*n);
    for i=1:n
        largeG(6*i-5:6*i,6*i-5:6*i) = robot.G(:,:,i);
    end 
    
    largeL = zeros(6*n,6*n);
    
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
            
        largeL(6*i-5:6*i,6*i-5:6*i) = eye(6,6);
        
        for j = 1:size(robot.tree{i}.children,1)
            stack_local = CStack();
            stack_local.push(robot.tree{i}.children(j));
            
            while(~stack_local.isempty)
                k = stack_local.pop;

                for child = 1:size(robot.tree{k}.children,1)
                    stack_local.push(robot.tree{k}.children(child));
                end
                
                largeL(6*i-5:6*i,6*k-5:6*k) = Ad_T(:,:,robot.tree{i}.children(j))'*largeL(6*robot.tree{i}.children(j)-5:6*robot.tree{i}.children(j),6*k-5:6*k);                
            end       
        end
    end
    
    largeL = largeL';

    M = largeA'*largeL'*largeG*largeL*largeA;
end