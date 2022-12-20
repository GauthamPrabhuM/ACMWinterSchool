clear
close all

rng(1);

prob = "burgers_1D";

load(prob + ".mat");    % Loads f, Z0 and xr variables from burgers1D_verification.m

steps = 500;

% dt = 1e-3; % Burgers
dt = 1e-1; % SWE
tspan = [0 dt];

atol = 1e-6;
rtol = 1e-6;
[t, u] = ode45(f, 0:tspan(2):steps*tspan(2), Z0', odeset('RelTol', rtol, 'AbsTol', atol));

%Plot this
figure(1);
plot(t,u);
xlabel("t");
ylabel("u");

figure(2);
plot(xr,u);
xlabel("x");
ylabel("u");
% legend (if needed)