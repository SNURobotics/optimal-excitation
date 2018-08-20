%% Caculate Derivatives of the Fourier Trajectory w.r.t. params
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]       [Description]                      [Size]
%  params       spline coefficients                2m*n    (parameter num * joints num)
%  w            base frequency (sin wt)            1*1
%  t            times to sample                    1*k     (1*desired number of frames)
%  q0           (optional) pos offset              n*1

%% Outputs
% [Name]       [Description]                                 [Size]
%  dq           spline pos derivative vector                  n*2m*n*k
%  dqdot        (optional) spline vel derivative vector       n*2m*n*k
%  dqddot       (optional) spline acc derivative vector       n*2m*n*k

%% Implementation
function [dq, dqdot, dqddot] = getFourierDerivative(params, w, t)
    n = size(params,2);
    m = floor(size(params,1)/2); % k = 1,2,...m
    
    num_t = size(t,2);
    
    dq = zeros(n, 2*m, n, num_t);
    dqdot = zeros(n, 2*m, n, num_t);
    dqddot = zeros(n, 2*m, n, num_t);
    
    for i = 1:n
        for j = 1:num_t
            for k = 1:m
                dq(i,2*k-1,i,j)     =   sin(k*w*t(j));
                dq(i,2*k,i,j)       =   cos(k*w*t(j));
                dqdot(i,2*k-1,i,j)  =   k*w*cos(k*w*t(j));
                dqdot(i,2*k,i,j)    = - k*w*sin(k*w*t(j));
                dqddot(i,2*k-1,i,j) = - k*k*w*w*sin(k*w*t(j));                
                dqddot(i,2*k,i,j)   = - k*k*w*w*cos(k*w*t(j));                
            end
        end
    end
end