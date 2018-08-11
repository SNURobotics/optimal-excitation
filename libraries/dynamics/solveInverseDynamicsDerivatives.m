%% Recursive Derivative Newton-Euler Inverse Dynamics Solver
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]       [Description]                      [Size]
%  A            i-th body screws from i-th frame   6*n    n = num of joints
%  M            initial relative frames M_{i,i-1}  4*4*n
%  q            joint angles                       n*1
%  qdot         joint vleocities                   n*1
%  qddot        joint accelerations                n*1
%  G            link inertial matrices             6*6*n
%  V            spatial velocities                 6*n
%  Vdot         spatial accelerations              6*n
%  dq           dq/dp                              n*m,   m = num of parameters
%  dqdot        dqdot/dp                           n*m
%  dqddot       dqddot/dp                          n*m
%  Vdot_0       (optional) base acceleration       6*1

%% Outputs
% [Name]  [Description]                                    [Size]
%  dtau    derivative of joint torques                      n*m
%  dV      (optional) derivative of spatial velocities      6*n*m
%  dVdot   (optional) derivative of spatial accelerations   6*n*m
%  dF      (optional) derivative of wrenches                6*n*m

%% Examples
% dtau = solveInverseDynamicsDerivatives(A,M,q,qdot,qddot,G,V,Vdot,dq,dqdot,dqddot)
% [dtau, dV, dVdot, dF] = solveInverseDynamicsDerivatives(A,M,q,qdot,qddot,G,V,Vdot,dq,dqdot,dqddot)

%% Implementation
function [dtau, varargout] = solveInverseDynamicsDerivatives(A,M,q,qdot,qddot,G,V,Vdot,dq,dqdot,dqddot,varargin)
    %% Initialization
    n      = size(q,1);         % number of joints
    m      = size(dq,2);        % number of parameters
    dV     = zeros(6,n,m);
    dVdot  = zeros(6,n,m);
    dF     = zeros(6,n,m);
    dtau   = zeros(n,m);
    
    V_0     = zeros(6,1);       % base velocity
    Vdot_0  = zeros(6,1);       % base acceleration
    dV_0    = zeros(6,1);       % derivative of base velocity
    dVdot_0 = zeros(6,1);       % derivative of base acceleration

    if     nargin == 11         % no optional inputs
    elseif nargin == 12         % optional base acceleration
        Vdot_0 = varargin{1};        
    else
        error('Init Error: undefined number of inputs');
    end
    
    T    = zeros(4,4,n); % T_{i,i-1}
    Ad_T = zeros(6,6,n); % Ad_T_{i,i-1}
    
    %% Forward Recursion
    for i = 1:n
        T(:,:,i)    = exp_se3(-A(:,i)*q(i))*M(:,:,i);
        Ad_T(:,:,i) = large_Ad(T(:,:,i));
        if i == 1
            for p = 1:m
                dV(:,i,p)    = Ad_T(:,:,i)*dV_0    - small_ad(A(:,i))*Ad_T(:,:,i)*V_0*dq(i,p)    + A(:,i)*dqdot(i,p);
                dVdot(:,i,p) = Ad_T(:,:,i)*dVdot_0 - small_ad(A(:,i))*Ad_T(:,:,i)*Vdot_0*dq(i,p) + small_ad(dV(:,i,p))*A(:,i)*qdot(i) ...
                                + small_ad(V(:,i))*A(:,i)*dqdot(i,p) + A(:,i)*dqddot(i,p);
            end
        else
            for p = 1:m
                dV(:,i,p)    = Ad_T(:,:,i)*dV(:,i-1,p)    - small_ad(A(:,i))*Ad_T(:,:,i)*V(:,i-1)*dq(i,p)    + A(:,i)*dqdot(i,p);
                dVdot(:,i,p) = Ad_T(:,:,i)*dVdot(:,i-1,p) - small_ad(A(:,i))*Ad_T(:,:,i)*Vdot(:,i-1)*dq(i,p) + small_ad(dV(:,i,p))*A(:,i)*qdot(i) ...
                                + small_ad(V(:,i))*A(:,i)*dqdot(i,p) + A(:,i)*dqddot(i,p);
            end
        end
    end

    %% Backward Recursion
    for i = n:-1:1
        if i == n
            for p = 1:m
                dF(:,i,p) = G(:,:,i)*dVdot(:,i,p) - small_ad(dV(:,i,p))'*G(:,:,i)*V(:,i) - small_ad(V(:,i))'*G(:,:,i)*dV(:,i,p);
            end
        else
            for p = 1:m
                dF(:,i,p) = Ad_T(:,:,i+1)'* dF(:,i+1,p) - Ad_T(:,:,i+1)'*small_ad(A(:,i+1))'*F(:,i+1)*dq(i+1,p) + G(:,:,i)*dVdot(:,i,p) ...
                            - small_ad(dV(:,i,p))'*G(:,:,i)*V(:,i) - small_ad(V(:,i))'*G(:,:,i)*dV(:,i,p);
            end
        end
        for p = 1:m
            dtau(i,p) = A(:,i)'*dF(:,i,p);
        end
    end
    
    varargout{1} = dV;
    varargout{2} = dVdot;
    varargout{3} = dF;    
end