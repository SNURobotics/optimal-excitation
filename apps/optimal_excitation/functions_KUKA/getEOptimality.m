%% Objective Function and its Derivative: E-optimality
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
%  f       1/lambda_min(C)                             1*1
%  grad    derivative of condition number of C(p)      m*n

%% Implementation
function [f, varargout] = getEOptimality(p, robot, trajectory, sigma_inv, varargin)
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
    f = eigenvalues(index_opt);
    if nargout > 1
        m = size(p,1);
        n = size(p,2);

        % grad cond(C)
        grad = zeros(m,n);
        Qt = Q';
        for i = 1:m
            for j =1:n
                grad(i,j) = - (Qt(index_opt, :) * gradC(:,:,i,j) * Q(:, index_opt))*(eigenvalues(index_opt)^2);
            end
        end
        varargout{1} = grad;
    end
end