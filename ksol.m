function K1 = ksol(n)
%KSOL Solves kernel equations, returns the control gain
%   n:  number of the u components in the n+1 system

y = linspace(1/n,1,n); % grid for y
xn = 257; x = linspace(0,1,xn); % grid for x
% continuum paramters (with constant lambda = mu = 1)
sigy = @(x,y,h) x.*(x+1).*(y-0.5).*x.^2.*(h-0.5);
Wy = @(x,y) x.*(x+1).*(y-0.5).*exp(x);
thy = @(x,y) -70*exp(x*35/pi^2).*y.*(y-1);
qy = @(y) cos(2*pi*y);
% precompute the function values into lookup tables (fewer function calls)
q = qy(y);
S = zeros(n,n);
WM = zeros(n,xn);
THM = WM;
S1 = zeros(n,n,xn);
S2 = S1;
for k=1:xn
  for s=1:n
    S(s,:) = sigy(x(k),y(1:n),y(s)); % S
  end
  % precompute thy and Wy
  THM(:,k) = thy(x(k),y(1:n));
  WM(:,k) = Wy(x(k),y(1:n));
  % precompute LU-decomposition of (2+S)
  [S1(:,:,k), S2(:,:,k)] = lu(2*eye(n)-S/(n*xn));
end
kd = zeros(1,1,n);
% kernels (k^i)_{i=1}^n and k^{n+1}
K = zeros(xn,xn,n); % 3D matrix, components x, xi, i
Kb = ones(xn,xn); % 2D matrix, compoents x, xi (initial guess)
% boundary condition for K on xi = x
for k = 1:xn
  for kn = 1:n
    K(k,k,kn) = -thy(x(k),y(kn))./2;
  end
end

% solve iteratively until changes are smaller than the chosen tolerance
iter = 400; % max number of iterations before termination
tol = 1e-8; % chosen tolerance to check solution convergemnce
dn = zeros(1,iter); % keeps track of changes in successive approximations
% k=x, m=xi
for mm = 1:iter
  % solve K (k^i for i = 1,2,...,n) based on Kb
  ok = K;
  okb = Kb;
  for k = 2:xn
    for m = 1:k-1
      % compute everything in 2D, then tranform to 3D
      % crude finite-difference approximation of the kernel equations
      kv = K(k-1,m,:) + K(k,m+1,:);
      kv = kv(:);
      kd(1,1,:) = S2(:,:,k)\(S1(:,:,k)\(kv + THM(:,m)*Kb(k,m)/xn));
      K(k,m,:) = kd;
    end
  end  
  % solve Kb (k^{n+1}) based on K
  Kb(:,1) = K(:,1,1)*q(1);
  for k = 2:n
    Kb(:,1) = Kb(:,1) + K(:,1,k)*q(k);
  end
  Kb(:,1) = Kb(:,1)/n;
  for k = 2:xn
    for m = 2:k-1
      % crude finite-difference approximation of the kernel equation
      kv = K(k,m,:);
      Kb(k,m) = 0.5*(Kb(k-1,m) + Kb(k,m-1) + sum(WM(:,m).*kv(:))/(n*xn));
    end
    % at the boundary xi=x, use m=k-1 as central point
    kv = K(k,k-1,:);
    Kb(k,k) = 2*Kb(k,k-1)-Kb(k-1,k-1) + sum(WM(:,k-1).*kv(:))/(n*xn);
  end
  % check changes, terminate if small, otherwise repeat
  dn(mm) = norm(K(:)-ok(:),inf) + norm(Kb(:)-okb(:),inf);
  if dn(mm) < tol
    % show terminating iteration
    disp(['Solution found on iteration ',num2str(mm)]) 
    break
  end
end
% evaluate control gains
Kn2 = zeros(n*(xn-1),1);
Knp = Kb(xn,1:xn-1)';
for k = 1:n
  Kn2((k-1)*(xn-1)+1:k*(xn-1)) = K(xn,2:xn,k)';
end
K1 = [Kn2; Knp];
end

