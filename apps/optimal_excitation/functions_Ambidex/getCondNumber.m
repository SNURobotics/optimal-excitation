%% Objective Function and its Derivative: Condition Number
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]      [Description]                                                [Size]
%  p           trajectory parameter to optimize                             m*n
%  robot       robot object (A, B, M)                                       struct
%  trajectory  object w/ trajectory parameters                              struct
%  sigma_inv   inverse of torque covariance matrix                          n*n
%  num_opt     (optional) specifed to get [num-opt]th largest eigenvalue    1*1  

%% Outputs
% [Name]  [Description]                               [Size]
%  f       cond(C)                                     1*1
%  grad    derivative of condition number of C(p)      m*n

%% Implementation
function [f, varargout] = getCondNumber(p, robot, trajectory, sigma_inv, varargin)
    num_opt = 1;    
    if nargin > 4
        num_opt = varargin{1};
    end
  
    if nargout == 1
        C = getObjectiveMatrixC(p, robot, trajectory, sigma_inv);
    else
        [C, gradC] = getObjectiveMatrixCwithGradient(p, robot, trajectory, sigma_inv);
    end
    
    % eigen decomposition
    [Q, D] = eig(C);
    D = inv(D);
    
    k = size(D,1);
    eigenvalues = zeros(k,1);
    for i=1:k
        eigenvalues(i) = D(i,i);    
    end
    
    [~, indicies] = sort(eigenvalues, 'descend');
    index_opt = indicies(num_opt,1);
    
    min_eig = min(eigenvalues);
    min_ind = find(eigenvalues == min_eig); min_ind = min_ind(1);

    
    f = eigenvalues(index_opt)/min_eig;
    if nargout > 1
        m = size(p,1);
        n = size(p,2);

        % grad cond(C)
        grad = zeros(m,n);
        Qt = Q';
        for i = 1:m
            for j =1:n
                grad(i,j) = (Qt(min_ind, :) * gradC(:,:,i,j) * Q(:, min_ind))*eigenvalues(index_opt) ...
                            - (Qt(index_opt, :) * gradC(:,:,i,j) * Q(:, index_opt))*eigenvalues(index_opt)^2/min_eig;
            end
        end
        varargout{1} = grad;
    end
end