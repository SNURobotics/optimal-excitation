%% Recursive Regressor Derivatives Calculator
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]       [Description]                        [Size]
%  A            i-th body screws from i-th frame     6*n    n = num of joints
%  M            initial relative frames M_{i,i-1}    4*4*n
%  q            joint angles                         n*1
%  V            spatial velocities                   6*n
%  Vdot         spatial accelerations                6*n
%  dq           spatial velocities                   6*n*m  m = num of parameters
%  dV           spatial velocities                   6*n*m  m = num of parameters
%  dVdot        spatial accelerations                6*n*m
%  W            (optional) sub-regressor matrix      6n*10n

%% Outputs
% [Name]  [Description]                                [Size]
%  dY      regressor derivative matrix                  n*10n*m
%  dW      (optional) sub-regressor derivative matrix   6n*10n*m                      

%% Implementation
function [dY, varargout] = getRegressorDerivativesRecursive(A,M,q,V,Vdot,dq,dV,dVdot,varargin)
    %% Initialization
    n      = size(q,1);         % number of joints
    m      = size(dq,2);        % number of parameters
    dW     = zeros(6*n, 10*n, m);
    
    if     nargin == 8          % no optional inputs
        [Y, W] = getRegressorRecursive(A,M,q,V,Vdot);
    elseif nargin == 9          % optional W
        W = varargin{1};
    else
        error('Init Error: undefined number of inputs');
    end
    
    diagA  = zeros(6*n,n);
    for i=1:n
        diagA(6*i-5:6*i,i) = A(:,i);
    end  
    
    T    = zeros(4,4,n); % T_{i,i-1}
    Ad_T = zeros(6,6,n); % Ad_T_{i,i-1}
    
    %% Recursion
    for i = n:-1:1
        T(:,:,i)    = exp_se3(-A(:,i)*q(i))*M(:,:,i);
        Ad_T(:,:,i) = large_Ad(T(:,:,i));
        
        for p=1:m
            dW(i,i,p) = convertVelocityToRegressor(dVdot(:,i,p)) - small_ad(dV(:,i,p))'*convertVelocityToRegressor(V(:,i)) ...
                         - small_ad(V(:,i))'*convertVelocityToRegressor(dV(:,i,p));
            for k = i+1:n
                dW(i,k,p) = Ad_T(:,:,i+1)'*dW(i+1,k,p) - Ad_T(:,:,i+1)'*small_ad(A(:,i+1))'*W(i+1,k)*dq(i+1,p);
            end
        end
    end
    
    dY = diagA'*dW;
    varargout{1} = dW;
end