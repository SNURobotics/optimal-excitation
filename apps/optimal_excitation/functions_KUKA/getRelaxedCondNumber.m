%% Objective Function and its Derivative: Relaxed Condition Number
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]      [Description]                                                [Size]
%  p           trajectory parameter to optimize                             m*n
%  robot       robot object (A, B, M)                                       struct
%  trajectory  object w/ trajectory parameters                              struct
%  sigma_inv   inverse of torque covariance matrix                          n*n

%% Outputs
% [Name]  [Description]                               [Size]
%  f       condition number of C(p)  10*1              1*1
%  grad    derivative of condition number of C(p)      m*n

%% Implementation
function [f, varargout] = getRelaxedCondNumber(p, robot, trajectory, sigma_inv)
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
    end
    
    % eigen decomposition
    tr_C = trace(C);
    C_inv = pinv(C); C_inv = (C_inv + C_inv')/2;
    tr_C_inv = trace(C_inv);
    f = tr_C*tr_C_inv / (size(C,1)^2);
    
    if nargout > 1
        m = size(p,1);
        n = size(p,2);

        % grad cond(C)
        grad = zeros(m,n);
        for i = 1:m
            for j =1:n
                grad(i,j) = (trace(gradC(:,:,i,j))*tr_C_inv - tr_C*trace(C_inv*gradC(:,:,i,j)*C_inv))/(size(C,1)^2);
            end
        end
        varargout{1} = grad;
    end
    eig_C = eig(C);
    cond = max(eig_C)/min(eig_C)
end