%% Trajectory Optimization for Parameter Identification
close all
clear
clc

%% Initialization
n           = 7;    % num of joints
m           = 10;   % num of spline coefs per joints
basis_order = 4;    % cubic spline
total_time  = 30;   % trajectory horizon

parameters = rand(m,n);

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fmincon %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = optimoptions(@fmincon,'Algorithm', 'sqp', 'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', false, 'TolCon',1e-5,'TolX',1e-10,'MaxFunEvals', ...
    1000000,'MaxIter',10000,'Display','iter','Hessian','bfgs'); %'fin-diff-grads'
% options = optimoptions(@fmincon,'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, 'Algorithm','sqp','TolCon',1e-15,'TolX',1e-15,'MaxFunEvals',1000000,'MaxIter',10000,'Display','iter');
% options = optimoptions(@fmincon, 'Algorithm','sqp','TolCon',1e-15,'TolX',1e-15,'MaxFunEvals',1000000,'MaxIter',10000,'Display','iter');
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
tic
[x,fval,exitflag,output,lam_costate] = fmincon(@getCondNumber, parameters, [], [], [], [], [], [], [], options);
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
