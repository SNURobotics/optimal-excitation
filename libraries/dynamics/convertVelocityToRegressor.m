%% convert velocity V to regressor y(V), where GV = y(V)phi
function V_regressor_mat = convertVelocityToRegressor(V)
V_regressor_mat = zeros(6,10);
w = V(1:3); v = V(4:6);
V_regressor_mat(4:6,1) = v;
V_regressor_mat(:,2:4) = [-skew(v);skew(w)];
V_regressor_mat(1:3,5:end) = [w(1),  0  ,  0   ,  w(2),   0  , w(3);
                               0  , w(2),  0   ,  w(1),  w(3),   0 ;
                               0  ,  0  , w(3) ,   0  ,  w(2), w(1) ];
end