%% Objective Function and its Derivative: Condition Number
% 2018 Bryan Dongik Lee

%% Inputs

%% Outputs

%% Implementation
function [f, grad] = getCondNumber(p)
    [C, gradC] = getObjectiveMatrixC(p);
    
    % eigen decomposition
    [Q, D] = eig(C);    
    k = size(D,1);
    eigs = zeros(k,1);
    for i=1:k
        eigs(i) = abs(D(i,i));    
    end
    min_eig = min(eigs);
    max_eig = max(eigs);
    min_ind = find(eigs == min_eig);
    max_ind = find(eigs == max_eig);

    % condition number cond(C)
    f = max_eig/min_eig;

    m = size(p,1);
    n = size(p,2);

    % grad cond(C)
    grad = zeros(m,n);
    Qt = Q';
    for i = 1:m
        for j =1:n
            grad(i,j) = (Qt(max_ind, :) * gradC(:,:,i) * Q(:, max_ind))/min_eig ...
                        - (Qt(min_ind, :) * gradC(:,:,i) * Q(:, min_ind))*max_eig/(min_eig^2);
        end
    end
end