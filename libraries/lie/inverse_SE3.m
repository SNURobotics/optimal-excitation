%% Analytic Inverse of SE(3)
function T_inv = inverse_SE3(T)

R = T(1:3,1:3);
p = T(1:3,4);
T_inv = [R.' -R.'*p; zeros(1,3) 1];