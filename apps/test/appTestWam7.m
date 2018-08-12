close all
clear
clc

%%
robot = wam7robot

n = robot.nDOF;

A = zeros(6,n);
for i=1:n
    A(:,i) = robot.link(i).screw;
end

M = zeros(4,4,n);
for i=1:n
    M(:,:,i) = inverse_SE3(robot.link(i).M);
end
