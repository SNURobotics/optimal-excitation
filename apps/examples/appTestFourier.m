close all
clear
clc

%%
params = [0.4 1; 3, 3]';
w = pi;
t = linspace(0,10,1000);
q0 = [10, 50]';

[q, qdot, qddot] = makeFourier(params, w, t, q0);

figure()
subplot(3,1,1);
plot(t,q(2,:))

subplot(3,1,2);
plot(t,qdot(2,:))

subplot(3,1,3);
plot(t,qddot(2,:))