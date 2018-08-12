%% convert velocity V to regressor y(V), where GV = y(V)phi
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]  [Description]                                     [Size]
%  V       spatial velocity                                  6*1

%% Outputs
% [Name]  [Description]                                     [Size]
%  y_V     linear mapping of velocity to regressor, y(V)     6*10

%% Implementation
function y_V = convertVelocityToRegressor(V)

    if(size(V) == [6, 1])
        y_V = zeros(6,10);
        w = V(1:3); v = V(4:6);
        y_V(4:6,1) = v;
        y_V(:,2:4) = [-skew(v);skew(w)];
        y_V(1:3,5:end) = [w(1),  0  ,  0   ,  w(2),   0  , w(3);
                                       0  , w(2),  0   ,  w(1),  w(3),   0 ;
                                       0  ,  0  , w(3) ,   0  ,  w(2), w(1) ];
    else
        error('convert error: wrong input size');
    end
end