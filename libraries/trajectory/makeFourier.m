%% Fourier Series Trajectory Generator a_k*sin(kwt) + b_k*cos(kwt)
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]       [Description]                      [Size]
%  params       spline coefficients                2m*n    (parameter num * joints num)
%  w            base frequency (sin wt)            1*1
%  t            times to sample                    1*k     (1*desired number of frames)
%  q0           (optional) pos offset              n*1

%% Outputs
% [Name]       [Description]                      [Size]
%  q            spline pos vector                  n*k
%  qdot         (optional) spline vel vector       n*k
%  qddot        (optional) spline acc vector       n*k

%% Implementation
function [q, qdot, qddot] = makeFourier(params, w, t, varargin)
    n = size(params,2);
    m = floor(size(params,1)/2); % k = 1,2,...m
    
    num_t = size(t,2);
    
    q = zeros(n, num_t);        
    qdot = zeros(n, num_t);
    qddot = zeros(n, num_t);
    
    if nargin > 3
        q = q + varargin{1};
    end
    
    for i = 1:n
        for j = 1:num_t
            for k = 1:m
                q(i,j) = q(i,j) + params(2*k-1,i)*sin(k*w*t(j)) + params(2*k,i)*cos(k*w*t(j));
                qdot(i,j) = qdot(i,j) + k*w*params(2*k-1,i)*cos(k*w*t(j)) - k*w*params(2*k,i)*sin(k*w*t(j));
                qddot(i,j) = qddot(i,j) - k*k*w*w*params(2*k-1,i)*sin(k*w*t(j)) - k*k*w*w*params(2*k,i)*cos(k*w*t(j));                
            end
        end
    end
end