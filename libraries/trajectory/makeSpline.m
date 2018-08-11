%% B-Spline Generator
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]       [Description]                      [Size]
%  params       spline coefficients                m*n     (parameter num * joints num)
%  order        order of spline basis functions    1*1
%  time         spline total horizon               1*1

%% Outputs
% [Name]       [Description]                      [Size]
%  sp           spline pos                         struct n
%  spdot        (optional) spline vel              struct n
%  spddot       (optional) spline acc              struct n

%% Implementation
function [sp, varargout] = makeSpline(params, order, time)
    n = size(params,2);
    m = size(params,1);
    
    num_knots = m + order;
    knots = zeros(1,num_knots);
    knots(1,num_knots-order+2:num_knots) = time;
    knots(1,order:num_knots-order+1) = linspace(0,time, num_knots-2*order+2);
    
    for i = 1:n
        sp(i) = spmak(knots, params(:,i)');
        spdot(i) = fnder(sp(i));
        spddot(i) = fnder(spdot(i));
    end
    
    varargout{1} = spdot;
    varargout{1} = spddot;    
end