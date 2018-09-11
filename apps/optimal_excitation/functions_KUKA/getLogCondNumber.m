%% Objective Function and its Derivative: Log Condition Number
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]      [Description]                                                [Size]
%  p           trajectory parameter to optimize                             m*n
%  robot       robot object (A, B, M)                                       struct
%  trajectory  object w/ trajectory parameters                              struct
%  sigma_inv   inverse of torque covariance matrix                          n*n

%% Outputs
% [Name]  [Description]                                      [Size]
%  f       log og condition number of C(p)  10*1              1*1
%  grad    derivative of log of condition number of C(p)      m*n

%% Implementation
function [f, varargout] = getLogCondNumber(p, robot, trajectory, sigma_inv)
    B_metric_inv_Phi_Bt_half = (robot.B_metric_inv_Phi_Bt)^(0.5);
    B_metric_inv_Phi_Bt_half = (B_metric_inv_Phi_Bt_half+B_metric_inv_Phi_Bt_half')/2;
    B_metric_inv_Phi_Bt_half_inv = pinv(B_metric_inv_Phi_Bt_half);
    B_metric_inv_Phi_Bt_half_inv = (B_metric_inv_Phi_Bt_half_inv+B_metric_inv_Phi_Bt_half_inv')/2;
    if nargout == 1
        C = getObjectiveMatrixC(p, robot, trajectory, sigma_inv);
        C = B_metric_inv_Phi_Bt_half*C*B_metric_inv_Phi_Bt_half_inv;
        C = (C+C')/2;
    else
        [C, gradC] = getObjectiveMatrixCwithGradient(p, robot, trajectory, sigma_inv);
        C = B_metric_inv_Phi_Bt_half*C*B_metric_inv_Phi_Bt_half_inv;
        C = (C+C')/2;
        for i = 1 : size(gradC,3)
            for j =1 : size(gradC,4)
                gradC(:,:,i,j) = B_metric_inv_Phi_Bt_half * gradC(:,:,i,j) * B_metric_inv_Phi_Bt_half_inv;
                gradC(:,:,i,j) = (gradC(:,:,i,j) + gradC(:,:,i,j)')/2;
            end
        end
        gradC = gradC_;
    end
    
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
    if nargout > 1
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
        varargout{1} = grad;
    end
end