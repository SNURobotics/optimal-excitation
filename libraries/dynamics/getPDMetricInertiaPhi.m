%% Compute Metric G(Phi) of Single Inertia Matrix Phi
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]  [Description]                                       [Size]
%  Phi     link inertia matrix for linear form dynamics        10*1

%% Outputs
% [Name]  [Description]                                       [Size]
%  metric  metric of inertia matrix phi                        10*10

%% Implementation
function metric = getPDMetricInertiaPhi(Phi)
    S = convertInertiaGToS(convertInertiaPhiToG(Phi));
    metric = zeros(10,10);
%     E = zeros(4,4,10);
%     E(4,4,1) = 1;
%     E(1,4,2) = 1; E(4,1,2) = 1;
%     E(2,4,3) = 1; E(4,2,3) = 1;
%     E(3,4,4) = 1; E(4,3,4) = 1;
%     E(1,1,5) = 1;
%     E(2,2,6) = 1;
%     E(3,3,7) = 1;
%     E(1,2,8) = 1; E(2,1,8) = 1;
%     E(2,3,9) = 1; E(3,2,9) = 1;
%     E(3,1,10) = 1; E(1,3,10) = 1;
    
    for j =1 : 10
        for k = 1: 10
            dphi_j = zeros(10,1); dphi_j(j) = 1;
            dphi_k = zeros(10,1); dphi_k(k) = 1;
            dP_j = convertInertiaGToS(convertInertiaPhiToG(dphi_j));
            dP_k = convertInertiaGToS(convertInertiaPhiToG(dphi_k));
            metric(j,k) = trace(pinv(S) * dP_j * pinv(S) * dP_k);
        end
    end
end