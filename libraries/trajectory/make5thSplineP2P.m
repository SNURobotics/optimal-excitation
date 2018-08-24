%% Point-to-Point 5-th Order Spline Generator
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]       [Description]                      [Size]
%  qi           initial point                      n*1
%  qidot        initial velocity                   n*1
%  qiddot       initial acceleration               n*1
%  qf           final velocity                     n*1
%  qfdot        final acceleration                 n*1
%  qfddot       final point                        n*1
%  horizon      total time horizon                 1*1
%  t            times to sample                    1*k     (1*desired number of frames)

%% Outputs
% [Name]       [Description]                      [Size]
%  q           spline pos                          n*k
%  qdot        (optional) spline vel               n*k
%  qddot       (optional) spline acc               n*k

%% Implementation
function [q, qdot, qddot] = make5thSplineP2P(qi, qidot, qiddot, qf, qfdot, qfddot, horizon, t)
    k = size(t,2);
    coef =[qi,...
        qidot,...
        qiddot/2,...
        -(20*qi - 20*qf + 8*qfdot*horizon + 12*qidot*horizon - qfddot*horizon^2 + 3*qiddot*horizon^2)/(2*horizon^3), ...
        (30*qi - 30*qf + 14*qfdot*horizon + 16*qidot*horizon - 2*qfddot*horizon^2 + 3*qiddot*horizon^2)/(2*horizon^4), ...
        -(12*qi - 12*qf + 6*qfdot*horizon + 6*qidot*horizon - qfddot*horizon^2 + qiddot*horizon^2)/(2*horizon^5)];
    time = [ones(1,k); t; t.^2; t.^3; t.^4; t.^5];
    
    q = coef*time;
    
    if nargout > 1
        timedot = [zeros(1,k); ones(1,k); 2*t; 3*t.^2; 4*t.^3; 5*t.^4];
        qdot = coef*timedot;
        timeddot = [zeros(1,k); zeros(1,k); 2*ones(1,k); 6*t; 12*t.^2; 20*t.^3];
        qddot = coef*timeddot;
    end
end