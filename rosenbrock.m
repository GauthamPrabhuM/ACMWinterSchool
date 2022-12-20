function f = rosenbrock(x,y)
% Introduced in 1960s
c = 1;
f = (1 - x).^2 + c*(y - x.^2).^2;
end

