n = 100;
A = hilb(n);
x = ones(n,1);
b = A*x;
x_solve = A\b;
error = norm(x_solve -x)