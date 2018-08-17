%% Objective Function and its Derivative: Condition Number
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]      [Description]                                                [Size]
%  p           trajectory parameter to optimize                             m*n
%  robot       robot object (A, B, M)                                       struct
%  trajectory  object w/ trajectory parameters(order, horizon, num_sample)  struct
%  sigma_inv   inverse of torque covariance matrix                          n*n

%% Outputs
% [Name]  [Description]                               [Size]
%  f       condition number of C(p)  10*1              1*1
%  grad    derivative of condition number of C(p)      m*n

%% Implementation
function [f, grad] = getLogCondNumber(p, robot, trajectory, sigma_inv)
    [C, gradC] = getObjectiveMatrixC(p, robot, trajectory, sigma_inv);
    
    % eigen decomposition
    [Q, D] = eig(C);    
    k = size(D,1);
    eigs = zeros(k,1);
    for i=1:k
        eigs(i) = abs(D(i,i));    
    end
    min_eig = min(eigs);
    max_eig = max(eigs);
    min_ind = find(eigs == min_eig); min_ind = min_ind(1);
    max_ind = find(eigs == max_eig); max_ind = max_ind(1);

    % condition number cond(C)
    f = log(max_eig/min_eig);

    m = size(p,1);
    n = size(p,2);

    % grad cond(C)
    grad = zeros(m,n);
    Qt = Q';
    for i = 1:m
        for j =1:n
            grad(i,j) = (Qt(max_ind, :) * gradC(:,:,i,j) * Q(:, max_ind))/max_eig ...
                        - (Qt(min_ind, :) * gradC(:,:,i,j) * Q(:, min_ind))/min_eig;
        end
    end
end