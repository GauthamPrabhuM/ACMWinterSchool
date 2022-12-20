function yp = lotka_eqns(t,y)
alpha = 1;
beta = 0.01;
gamma = 0.02;
delta = 1;
yp = diag([alpha - beta*y(2), -delta + gamma*y(1)])*y;
end
