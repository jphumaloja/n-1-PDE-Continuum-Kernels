function [usol, Usol] = n1sim(n,T,tg,Ke)
%N1SIM Simulates the n+1 system for t=[0,T]
% return solution with timegrid linspace(0,T,tg)
%   n:  number of the u components in the n+1 system
%   T:  simulation end time
%   tg: timegrid (number of points) to evaluate the solution on
%   Ke: user-provided control gain

% continuum paramters
lamy = @(x,y) ones(size(x));
muy = @(x) ones(size(x));
sigy = @(x,y,h) x.*(x+1).*(y-0.5).*x.^2.*(h-0.5);
Wy = @(x,y) x.*(x+1).*(y-0.5).*exp(x);
thy = @(x,y) -70*exp(x*35/pi^2).*y.*(y-1);
qy = @(y) cos(2*pi*y);

% parameters for computations
y = linspace(1/n,1,n); % grid for y
xg = 256; % grid for x (minus one end point determined by boundary conditions)
tmp = linspace(0,1,xg+1); % auxiliary grid
x0 = tmp(1:xg); x1 = tmp(2:xg+1); % grids for u and v

% n+1 system parameters
lam = @(i,x) lamy(x,y(i));
mu = @(x) muy(x);
sig = @(i,j,x) sigy(x,y(i),y(j));
W = @(i,x) Wy(x,y(i));
th = @(j,x) thy(x,y(j));
q = qy(y);

% finite difference approximation \dot z = Az + BU, where
A = zeros((n+1)*xg); % system operator initialization)
B = [zeros((n+1)*xg-1, 1); mu(1)*xg]; % control operator; v(1) = U

% auxiliary: D is backward difference, -D' is forward difference
D = eye(xg) - diag(ones(1,xg-1), -1);

% fill A
for k = 1:n
  Ik = (k-1)*xg+1:k*xg; % index set
  A(Ik,Ik) = -xg*lam(k,x1).*D; % tranport term
  A(Ik(1),n*xg+1) = xg*lam(k,0)*q(k); % boundary condition
  for l = 1:n
    Il = (l-1)*xg+1:l*xg; % index set
    A(Ik,Il) = A(Ik,Il) + diag(sig(k,l,x1))/n; % sigma terms (scale 1/n)
  end
  Il = n*xg+1:(n+1)*xg; % index set
  A(Ik,Il) = A(Ik,Il) + diag(W(k,x0)); % W term
end
k=n+1;
Ik = (k-1)*xg+1:k*xg; % index set
A(Ik,Ik) = -xg*mu(x0).*D'; % transport term
for l = 1:n
  Il = (l-1)*xg+1:l*xg; % index set
  A(Ik,Il) = A(Ik,Il) + diag(th(l,x1))/n; % theta terms (scale 1/n)
end

if nargin < 4 % use continuum kernels, if no user-provided gain
% kernel solutions and control law
kb = @(x,xi) 35/(2*pi^2);
ky = @(x,xi,y) 35*y.*(y-1).*exp(2*xi.*kb(x,xi));
% vectorize ky, discretize to n components, kb as is (the n+1 part)
Kn = zeros(xg*n, 1);
for k = 1:n
  Kn((k-1)*xg+1:k*xg) = ky(1, x1, y(k))';
end
% the control law (ode-solver compatible form)
U = @(z) (trapz(Kn.*z(1:n*xg)/n) + trapz(kb(1,x0).*z(n*xg+1:(n+1)*xg)))/xg;
else % use user-provided control gain
Ke2 = Ke(xg*n+1:end);
Ke1 = Ke(1:xg*n); 
U = @(z) (trapz(Ke1.*z(1:n*xg)/n) + trapz(Ke2.*z(n*xg+1:(n+1)*xg)))/xg;
end

z0 = [kron(q', ones(xg,1)); ones(xg,1)]; % initial condition
% simulate for [0,T] with initial condition z0, return result
opts = odeset('Jacobian', A); % should expedite computations
% solve; ode45 (ok fast), ode23 (faster, accuracy?), others much slower
sol = ode45(@(t,z) A*z + B*U(z), [0, T], z0, opts);
TT = linspace(0,T,tg); % time grid (for plotting)
usol = deval(sol, TT); % evaluate solution at tg
Usol = zeros(1,tg); % evaluate input at tg
for k=1:tg % compute and return inputs
  Usol(k) = U(usol(:,k));
end
end

