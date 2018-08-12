%% Closed-Form Forward Dynamics Solver
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]       [Description]                      [Size]
%  A            i-th body screws from i-th frame   6*n
%  M            initial relative frames M_{i,i-1}  4*4*n
%  q            joint angles                       n*1
%  qdot         joint vleocities                   n*1
%  tau          joint torques                      n*1
%  G            link inertial matrices             6*6*n
%  Vdot_0       (optional1) base acceleration      6*1

%% Outputs
% [Name]       [Description]                      [Size]
%  qddot        joint accelerations                n*1

%% Implementation
function qddot = solveForwardDynamics(A,M,q,qdot,tau,G,varargin)
    %% Initialization
    n     = size(q,1);          % number of joints
    
    V_0    = zeros(6,1);        % base velocity
    Vdot_0 = zeros(6,1);        % base acceleration
    if     nargin == 6
    elseif nargin == 7
        Vdot_0  = varargin{1};  % optional base acceleration
    end
  
    T    = zeros(4,4,n);        % T_{i,i-1}
    Ad_T = zeros(6,6,n);        % Ad_T_{i,i-1}

    %% A
    diagA  = zeros(6*n,n);
    for i=1:n
        diagA(6*i-5:6*i,i) = A(:,i);
    end
    
    %% G
    diagG = zeros(6*n,6*n);
    for i=1:n
        diagG(6*i-5:6*i,6*i-5:6*i) = G(:,:,i);
    end
    
    %% ad_V
    V    = zeros(6,n);
    ad_V = zeros(6*n,6*n);
    for i = 1:n
        T(:,:,i)    = exp_se3(-A(:,i)*q(i))*M(:,:,i);
        Ad_T(:,:,i) = large_Ad(T(:,:,i));
        if i == 1
            V(:,i) = Ad_T(:,:,i)*V_0  + A(:,i)*qdot(i);
        else
            V(:,i) = Ad_T(:,:,i)*V(:,i-1)  + A(:,i)*qdot(i);
        end
        
        ad_V(6*i-5:6*i,6*i-5:6*i) = small_ad(V(:,i));
    end

    %% Vdot_base
    Vdot_base = zeros(6*n,1);
    Vdot_base(1:6,1) = Ad_T(:,:,1)* Vdot_0;
    
    %% ad_A_qdot
    ad_A_qdot = zeros(6*n,6*n);
    for i=1:n
        ad_A_qdot(6*i-5:6*i,6*i-5:6*i) = small_ad(A(:,i)*qdot(i));
    end
    
    %% W
    W = zeros(6*n,6*n);
    for i = 2:n
        W(6*i-5:6*i,6*i-11:6*i-6) = Ad_T(:,:,i);
    end
    
    %% L
    L = eye(6*n);
    for i=2:n
        for j=i-1:-1:1
            L(6*i-5:6*i,6*j-5:6*j) = Ad_T(:,:,i)*L(6*i-11:6*i-6,6*j-5:6*j);
        end
    end
    
    %% Closed-Form Dynamics
    M_q = diagA'*L'*diagG*L*diagA;
    c_q = -diagA'*L'*(diagG*L*ad_A_qdot*W + ad_V'*diagG)*L*diagA*qdot;
    g_q = diagA'*L'*diagG*L*Vdot_base;
    
    qddot = M_q\(tau - c_q - g_q);
end