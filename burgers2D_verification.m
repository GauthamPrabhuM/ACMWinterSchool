clear
close all

rng(3);

nx = 5;
x = linspace(0, 1, nx + 2).';
x = x(2:end-1);
dx = x(2) - x(1);

ny = nx;
y = linspace(0, 1, ny + 2);
y = y(2:end-1);
dy = y(2) - y(1);

[X, Y] = ndgrid(x, y);

n = nx*ny*2;
bval = 3;

Z0 = zeros(nx, ny, 2);
Z0(:, :, 1) = sin(pi*x).*sin(pi*y);
Z0(:, :, 2) = sin(pi*x).*sin(pi*y);
Z0 = Z0(:);

% binds = (Z0 == bval);
% Z0(binds) = 0;
% Z0 = Z0(:);

ka = 1.4;
kb = 0.1;

[Dx, Dy, Lxy] = diff(nx, ny, dx, dy);

res = @(Y) deal( Y(1:n/2, :), Y(1+n/2:n, :));
res2d = @(Y1) reshape(Y1, nx, ny);

f = @(t, Y) rhs(Y, res, Dx, Dy, Lxy, ka, kb);

save('burgers_2D.mat', 'f', 'Z0');
% J = @(t, Y) rhs_Y(Y, res, Dx, Dy, Lxy, ka, kb, n/2);
% Ja = @(t, Y) rhs_Ya(Y, res, Dx, Dy, Lxy, ka, kb, n/2);
% Jv = @(t, Y, vx) rhs_Yv(vx, Y, res, Dx, Dy, Lxy, ka, kb);
% Jav = @(t, Y, vx) rhs_Yav(vx, Y, res, Dx, Dy, Lxy, ka, kb);
% 
% [indsx, indsy] = ndgrid(x, y);
% indsx = indsx(:);
% indsy = indsy(:);
% B1 = 1e-1*exp(-sqrt((indsx - indsx.').^2 + (indsy - indsy.').^2)*2);
% B = [B1, 0.7*B1; 0.7*B1, B1];
% sqrtmB = sqrtm(B);
% Bi = inv(B);

% syms vv [n, 1]
% matlabFunction(B*vv, 'Vars', {vv}, 'Optimize', false, 'File', '+gen/b2Bvfunc');
% matlabFunction(Bi*vv, 'Vars', {vv}, 'Optimize', false, 'File', '+gen/b2Bivfunc');
% matlabFunction(sqrtmB*vv, 'Vars', {vv}, 'Optimize', false, 'File', '+gen/b2sqrtmBvfunc');
% Bv = @(v) gen.b1Bvfunc(v);
% Biv = @(v) gen.b1Bivfunc(v);
% sqrtmBv = @(v) gen.b1sqrtmBvfunc(v);

% Bv = @(v) B*v;
% Biv = @(v) Bi*v;
% sqrtmBv = @(v) sqrtmB*v;
% 
% xt = Z0;
% xf = sqrtmB*randn(n, 1) + xt;

% obsfac = 0.2;
% inds = randperm(length(nbinds));
% inds = inds(1:round(obsfac*n/2));
% inds = round(linspace(1, n/2, obsfac*n/2)).';
% obsvars = [inds; inds + n/2];
% 
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


% Z = Z0 + 10*randn(n, 1);
% jac = J(0, Z);
% jacT = Ja(0, Z);
% eyee = speye(n, n);
% 
% jvpi = sparse(n, n);
% javpi = sparse(n, n);
% fdjac = sparse(n, n);
% 
% epsi = 1e-6;
% for i = 1:n
%     fdjac(:, i) = (f(0, Z + epsi*eyee(:, i)) - f(0, Z))/epsi;
%     jvpi(:, i) = Jv(0, Z, eyee(:, i));
%     javpi(:, i) = Jav(0, Z, eyee(:, i));
% end
% 
% sqrt(sum((fdjac - jac).^2, 'all'))
% jac - jacT.'
% jac - jvpi
% jacT - javpi

% ax1 = surf(X, Y, reshape(Z0(1:n/2), nx, ny), 'EdgeColor', 'none');
% view(2)
% zlim([0 1])
% colorbar
% caxis([0 1])
% 
% dt = 1e-3;
% Z = Z0;
% for i = 1:1e3
%     
%     [~, Z] = DOPRI54(f, [0 dt], Z, 1e-3, 1e-3);
%     
%     ax1.ZData = reshape(Z(1:n/2), nx, ny);
%     title("Title = "+i*dt);
%     drawnow;
% 
% end

% 
function dY = rhs(Y, rf, Dx, Dy, Lxy, a, b) 

[u, v] = rf(Y);

u_t = - a*(u.*(Dx*u) + v.*(Dy*u)) + b*(Lxy*u);
v_t = - a*(u.*(Dx*v) + v.*(Dy*v)) + b*(Lxy*v);

dY = [u_t; v_t];

end
 
% 
% function dY = rhs_Y(Y, rf, Dx, Dy, Lxy, a, b, n) 
% 
% [u, v] = rf(Y);
% 
% u_x = Dx*u;
% u_y = Dy*u;
% 
% v_x = Dx*v;
% v_y = Dy*v;
% 
% uDxvDy = u.*Dx + v.*Dy;
% 
% u_tu = - a*(spdiags(u_x, 0, n, n) + uDxvDy) + b*(Lxy);
% u_tv = - a*spdiags(u_y, 0, n, n);
% 
% v_tu = - a*spdiags(v_x, 0, n, n);
% v_tv = - a*(spdiags(v_y, 0, n, n) + uDxvDy) + b*(Lxy);
% 
% dY = [u_tu, u_tv; v_tu, v_tv];
% 
% end
% 
% function dY = rhs_Ya(Y, rf, Dx, Dy, Lxy, a, b, n)
% 
% [u, v] = rf(Y);
% 
% u_x = Dx*u;
% u_y = Dy*u;
% 
% v_x = Dx*v;
% v_y = Dy*v;
% 
% uDxvDyT = (u.*Dx + v.*Dy).';
% 
% u_tu = - a*(spdiags(u_x, 0, n, n) + uDxvDyT) + b*(Lxy.');
% u_tv = - a*spdiags(u_y, 0, n, n);
% 
% v_tu = - a*spdiags(v_x, 0, n, n);
% v_tv = - a*(spdiags(v_y, 0, n, n) + uDxvDyT) + b*(Lxy.');
% 
% dY = [u_tu, v_tu; u_tv, v_tv;];
% 
% end
% 
% function dY = rhs_Yv(vx, Y, rf, Dx, Dy, Lxy, a, b) 
% 
% [u, v] = rf(Y);
% [vxu, vxv] = rf(vx);
% 
% u_x = Dx*u;
% u_y = Dy*u;
% 
% v_x = Dx*v;
% v_y = Dy*v;
% 
% u_tvx = - a*(u_x.*vxu + u.*(Dx*vxu) + v.*(Dy*vxu) + u_y.*vxv) + b*(Lxy*vxu);
% v_tvx = - a*(v_x.*vxu + v_y.*vxv + u.*(Dx*vxv) + v.*(Dy*vxv)) + b*(Lxy*vxv);
% 
% dY = [u_tvx; v_tvx];
% 
% end
% 
% function dY = rhs_Yav(vx, Y, rf, Dx, Dy, Lxy, a, b) 
% 
% [u, v] = rf(Y);
% [vxu, vxv] = rf(vx);
% 
% u_x = Dx*u;
% u_y = Dy*u;
% 
% v_x = Dx*v;
% v_y = Dy*v;
% 
% u_tvx = - a*(u_x.*vxu + Dx.'*(u.*vxu) + Dy.'*(v.*vxu) + v_x.*vxv) + b*(Lxy.'*vxu);
% v_tvx = - a*(u_y.*vxu + v_y.*vxv + Dx.'*(u.*vxv) + Dy.'*(v.*vxv)) + b*(Lxy.'*vxv);
% 
% dY = [u_tvx; v_tvx];
% 
% end
% 
% 
function [Dx, Dy, Lxy] = diff(nx, ny, dx, dy)

e1 = (0.5/dx)*ones(nx, 1);
Dx = spdiags([-e1, e1], [-1, 1], nx, nx);

e1 = (1/dx/dx)*ones(nx, 1);
Lx = spdiags([e1, -2*e1, e1], [-1, 0, 1], nx, nx);

e1 = (0.5/dy)*ones(ny, 1);
Dy = spdiags([e1, -e1], [-1, 1], ny, ny);

e1 = (1/dy/dy)*ones(ny, 1);
Ly = spdiags([e1, -2*e1, e1], [-1, 0, 1], ny, ny);

Dx = kron(speye(ny), Dx);
Dy = kron(Dy.', speye(nx));
Lxy = kron(Ly.', speye(nx)) + kron(speye(ny), Lx);

end
