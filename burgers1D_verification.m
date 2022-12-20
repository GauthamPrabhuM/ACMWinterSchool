clear
rng(3);

n = 100;
x = linspace(0, 1, n + 2).';    % Select 100 equidistant points between 0 & 1
dx = x(2) - x(1);

xr = x(2:end-1);
Z0 = [sin(pi*xr)];

ka = 1.4;
kb = 0.001;

[Dx, Lx] = diff(n, dx);

f = @(t, u) - ka*(u.*(Dx*u)) + kb*(Lx*u);

save('burgers_1D.mat','Z0', 'f', 'xr')
% f = @(t, u) - 0.5*randka*(Dx*(u.*u)) + kb*(Lx*u);
% J = @(t, u) - ka*(spdiags(Dx*u, 0, n, n) + spdiags(u, 0, n, n)*Dx) + kb*Lx;
% Jv = @(t, u, v) - ka*((Dx*u).*v + u.*(Dx*v)) + kb*(Lx*v);
% Jav = @(t, u, v) - ka*((Dx*u).*v + Dx.'*(u.*v)) + kb*(Lx.'*v);
% 
% inds = 1:n;
% B = 1e-1*exp(-abs(inds - inds.')/50);
% sqrtmB = sqrtm(B);
% Bi = inv(B);
% 
% randvec = randn(n, 1);
% 
% syms vv [n, 1]
% matlabFunction(B*vv, 'Vars', {vv}, 'Optimize', false, 'File', '+gen/b1Bvfunc');
% matlabFunction(Bi*vv, 'Vars', {vv}, 'Optimize', false, 'File', '+gen/b1Bivfunc');
% matlabFunction(sqrtmB*vv, 'Vars', {vv}, 'Optimize', false, 'File', '+gen/b1sqrtmBvfunc');
% 
% %%% This was for experimentation
% 
% Bv = @(v)B*v;
% Biv = @(v)B\v;
% sqrtmBv = @(v)sqrtmB*v
% 
% %Bv = @(v) gen.b1Bvfunc(v);
% %Biv = @(v) gen.b1Bivfunc(v);
% %sqrtmBv = @(v) gen.b1sqrtmBvfunc(v);
% 
% xt = Z0;
% xf = sqrtmBv(randvec) + xt;
% 
% obsvars = round(linspace(1, n, 20));
% nobs = length(obsvars);
% 
% H = eye(n);
% H = H(obsvars, :);
% 
% sigR = 1e-2;
% R = sigR*speye(nobs);
% sqrtR = sqrt(sigR)*speye(nobs);
% Ri = (1/sigR)*speye(nobs);
% 
% save('burgers1D.mat', 'f', 'J', 'Jv', 'Jav', 'n', 'obsvars', 'nobs', 'H', 'R', 'sqrtR', 'Ri', 'Bv', 'Biv', 'sqrtmBv', 'xt', 'xf')
% 
% % Z = Z0 + 10*randn(n, 1);
% % jac = J(0, Z);
% % jacT = jac.';
% % eyee = speye(n, n);
% % 
% % jvpi = sparse(n, n);
% % javpi = sparse(n, n);
% % fdjac = sparse(n, n);
% % 
% % epsi = 1e-6;
% % for i = 1:n
% %     fdjac(:, i) = (f(0, Z + epsi*eyee(:, i)) - f(0, Z))/epsi;
% %     jvpi(:, i) = Jv(0, Z, eyee(:, i));
% %     javpi(:, i) = Jav(0, Z, eyee(:, i));
% % end
% % 
% % sqrt(sum((fdjac - jac).^2, 'all'))
% % jac - jacT.'
% % jac - jvpi
% % jacT - javpi
% 
% % ax1 = plot(x, [0; Z0; 0]);
% % xlim([0 1]);
% % ylim([0 1]);
% % dt = 1e-3;
% % Z = Z0;
% % for i = 1:1000
% %     [~, Z] = DOPRI54(f, [0 dt], Z, 1e-3, 1e-3);
% %     ax1.YData(2:end-1) = Z;
% %     drawnow;
% % end
% % 
function [Dx, Lx] = diff(n, dx)
 
e1 = (0.5/dx)*ones(n, 1);
Dx = spdiags([-e1, e1], [-1, 1], n, n);
 
e1 = (1/dx/dx)*ones(n, 1);
Lx = spdiags([e1, -2*e1, e1], [-1, 0, 1], n, n);
 
end
