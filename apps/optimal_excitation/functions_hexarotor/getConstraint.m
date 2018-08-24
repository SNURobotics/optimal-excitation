%% Constraints and its Derivatives for Trajectory Optimization
% 2018 Bryan Dongik Lee

%% Inputs

%% Outputs

%% Implementation
function [c, ceq, gradc, gradceq] = getConstraint(p, trajectory, robot)
    ceq = [];

    p_max = [1e10 1e10 1e10 1e10 1e10 1e10
              3  3  3 4*pi 4*pi 4*pi
              1  1  1  2  2  2];
    p_min = -p_max;
    c = [p - p_max; p_min - p];
    
    if nargout > 2
        gradceq = [];
        gradc = [ones(size(p));-ones(size(p))];
    end
    
%     % input constraints
%     u_min = robot.u_min;
%     u_max = robot.u_max;
%     
%     % Spline trajectory parameters
%     num_sample       = trajectory.num_sample;
%     horizon          = trajectory.horizon;
%     trajectory_order = trajectory.order;
%     sample_time      = linspace(0,horizon,num_sample);
%     
%     % Fourier trajectory parameters
% %     num_sample       = trajectory.num_sample;
% %     horizon          = trajectory.horizon;
% %     base_frequency   = trajectory.base_frequency;
% %     sample_time      = linspace(0,horizon,num_sample);
% 
%     n = robot.dof;
%     G = zeros(6,6,n);                    % dummy value
%     S_inv = robot.S_inv;
%     Vdot_0 = [0;0;0;0;0;9.8];
% 
%     c = zeros(2*n,num_sample);
% 
%     % Spline trajectory generation with parameter p
%     [q, qdot, qddot]    = makeSpline(p, trajectory_order, horizon, sample_time);
% 
%     % Spline trajectory generation with parameter p
% %     [q, qdot, qddot]    = makeFourier(p, base_frequency, sample_time);
% 
%     if nargout == 2
%         for t = 1:num_sample
%             [tau, V, Vdot] = solveInverseDynamics(robot.A,robot.M,q(:,t),qdot(:,t),qddot(:,t),G,Vdot_0);
%             Y = S_inv*(convertVelocityToRegressor(Vdot(:,n)) - small_ad(V(:,n))'*convertVelocityToRegressor(V(:,n)));        
%             u = Y*robot.Phi;
%             c(:,t) = [u-u_max; u_min-u];
%         end
%     elseif naragout > 2
%         gradceq = [];
%         m = size(p,1);
%         
%         [dq, dqdot, dqddot] = getSplineDerivative(p, trajectory_order, horizon, sample_time);
%         dY = zeros(n,10,m,n);
% 
%         gradc = zeros(2*n,num_sample,m,n);
%         
%         for t = 1:num_sample
%             [tau, V, Vdot] = solveInverseDynamics(robot.A,robot.M,q(:,t),qdot(:,t),qddot(:,t),G,Vdot_0);
%             [dtau, dV, dVdot] = solveInverseDynamicsDerivatives(A,M,q(:,t),qdot(:,t),G,V,Vdot,dq(:,:,:,t),dqdot(:,:,:,t),dqddot(:,:,:,t),F,Vdot_0);
% 
%             Y = S_inv*(convertVelocityToRegressor(Vdot(:,n)) - small_ad(V(:,n))'*convertVelocityToRegressor(V(:,n)));
%             u = Y*robot.Phi;
%             c(:,t) = [u-u_max; u_min-u];
%             
%             for p = 1:m
%                 for l = 1:n              
%                     du = S_inv*(convertVelocityToRegressor(dVdot(:,n,p,l)) - small_ad(dV(:,n,p,l))'*convertVelocityToRegressor(V(:,n)) ...
%                                          - small_ad(V(:,n))'*convertVelocityToRegressor(dV(:,n,p,l)))*robot.Phi;
%                     gradc(:,t,p,l) = [du; -du];
%                 end
%             end            
%         end
%     end
end