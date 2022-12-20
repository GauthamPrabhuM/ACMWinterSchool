% Lotka Volterra Model
clear;
clc;
close all;
%% Writing Equations
% Firstly, you need the differential equations you are trying to solve
% (NOTE: See file lotka_eqns.m)

%% Simulate System
% After having defined your differential equations, you just need to decide
% on first, what numerical differential solver you would like to solve.

% We'll set the time interval: 0 < t < 15
t0 = 0;
tfinal = 15;

% Defining initial conditions: Let's say we have 20 units of population of
% both predators and prey
y0 = [20; 20];

% Solver: (Using ode23)
tic;
[t23,y23] = ode23(@lotka_eqns,[t0,tfinal],y0);
toc;

tic;
[t15s,y15s] = ode15s(@lotka_eqns,[t0,tfinal],y0);
toc;
%% Plotting (Of both numerical methods)
% % Plotting both populations with time
% figure(1)
% plot(t23,y23)
% hold on;
% plot(t15s,y15s)
% title('Predator/Prey Populations Over Time')
% xlabel('t')
% ylabel('Population')
% legend('Prey','Predators','Location','North')

% Plotting both populations against each other
figure(1)
plot(y23(:,1),y23(:,2))
hold on;
plot(y15s(:,1),y15s(:,2))
title('Phase Plane Plot')
xlabel('Prey Population')
ylabel('Predator Population')
legend('ode23','ode15s')