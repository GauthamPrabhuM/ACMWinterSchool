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
tfinal = 50;

% Defining initial conditions: Let's say we have 20 units of population of
% both predators and prey
y0 = [20; 20];

% Solver: (Using ode23)
[t,y] = ode23(@lotka_eqns,[t0,tfinal],y0);

%% Plotting
% Plotting both populations with time
figure(1)
plot(t,y)
title('Predator/Prey Populations Over Time')
xlabel('t')
ylabel('Population')
legend('Prey','Predators','Location','North')

%% Plotting both populations against each other
figure(2)
plot(y(:,1),y(:,2))
title('Phase Plane Plot')
xlabel('Prey Population')
ylabel('Predator Population')