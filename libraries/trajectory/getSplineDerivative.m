%% Calculate Derivative of B-Spline w.r.t. Parameters
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]       [Description]                      [Size]
%  params       spline coefficients                m*n     (parameter num * joints num)
%  order        order of spline basis functions    1*1
%  horizon      total time horizon                 1*1
%  t            (optional) user specified times    1*t     (1*desired number of frames)

%% Outputs
% [Name]       [Description]                                  [Size]
%  dsp           spline pos derivative functions               struct n*m*n
%  dspdot        (optional) spline vel derivative functions    struct n*m*n
%  dspddot       (optional) spline acc derivative functions    struct n*m*n
%  If t is specified, return real-valued vectors   n*m*n*t

%% Implementation
function [dsp, varargout] = getSplineDerivative(params, order, horizon, varargin)
    n = size(params,2);
    m = size(params,1);
    
    num_knots = m + order;
    knots = linspace(0, horizon, num_knots);
    
    for i = 1:n
        for j = 1:m
            p_i = zeros(m,1);
            p_i(j) = 1;
            dsp(i,j,i)     = spmak(knots, p_i');
            dspdot(i,j,i)  = fnder(dsp(i,j,i));
            dspddot(i,j,i) = fnder(dspdot(i,j,i)); 
        end
    end
    
    if     nargin == 3
    elseif nargin == 4   % if time specified
        t = varargin{1};
        num_t = size(t,2);
        dq = zeros(n, m, n, num_t);
        dqdot = zeros(n, m, n, num_t);
        dqddot = zeros(n, m, n, num_t);
        
        for i = 1:n
            for j = 1:m
                dq(i,j,i,:)     = fnval(dsp(i,j,i), t);
                dqdot(i,j,i,:)  = fnval(dspdot(i,j,i), t);
                dqddot(i,j,i,:) = fnval(dspddot(i,j,i), t);
            end
        end

        dsp = dq;
        dspdot = dqdot;
        dspddot = dqddot;
    end
    
    if nargout > 1
        varargout{1} = dspdot;
        varargout{2} = dspddot;
    end
end