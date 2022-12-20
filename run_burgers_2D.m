clear
close all

rng(1);

prob = "burgers_2D";

load(prob + ".mat");    % Loads f and Z0 from burgers2D_verification.m

steps = 100;

% dt = 1e-3; % Burgers
dt = 1e-2; % SWE
tspan = [0 dt];

atol = 1e-6;
rtol = 1e-6;
[t, u] = ode45(f, 0:tspan(2):steps*tspan(2), Z0', odeset('RelTol', rtol, 'AbsTol', atol));

%Plot this
figure(1);
plot(t,u);
xlabel("t");
ylabel("u");

% figure(2); (substitute for x in 2D case ?)
% plot(x,u);
% xlabel("x");
% ylabel("u");
% legend (if needed)