% Rosenbrock
clear;
clc;
close all;
    
%% Writing Function Equation
% Function has been defined in rosenbrock.m file

%% Optimize using MATLAB's Optimization toolbox
x = optimvar('x',1,2);
obj = rosenbrock(x(1),x(2));
prob = optimproblem('Objective',obj);

show(prob);

x0.x = [-5 5];
[sol,fval,exitflag,output] = solve(prob,x0);
%% Descent animation / Plotting
figure(1);
[X,Y] = meshgrid(-5:0.1:5,-5:0.1:5);
z = rosenbrock(X,Y);
surf(X,Y,z);
shading interp;
hold on;
plot3(sol.x(1),sol.x(2),rosenbrock(sol.x(1),sol.x(2)),'ro');
hold on;
plot3(x0.x(1),x0.x(2),rosenbrock(x0.x(1),x0.x(2)),'go');
xlabel('X axis');
ylabel('Y axis');
zlabel('Z axis');
