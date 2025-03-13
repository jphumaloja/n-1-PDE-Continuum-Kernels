%% Collect simulation data
nn = 2:6; % selected values of n
% initialize
cu = cell(size(nn));
cku = cu;
% generic simulation parameters
T = 10; % simulation end time
xg = 256; % x grid
tmp = linspace(0,1,xg+1); % auxiliary grid
x0 = tmp(1:xg); x1 = tmp(2:xg+1); % grids for plotting u and v
tg = 513; % time grid (for plotting the solution)
TT = linspace(0,T,tg); % time vector of solution
% simulate for fixed initial conditions (depending on n)
for k = 1:numel(nn)
  [usol, Usol] = n1sim(nn(k),T,tg);
  cu{k} = usol;
  cku{k} = Usol;
end
%% input plot (Fig. 2)
f1 = figure(1);
f1.Position = 40*[5 5 15 7];
for k = 1:numel(nn)
  lbl = ['$n=', num2str(nn(k)), '$'];
  plot(TT, cku{k},'-','linewidth', 2,'DisplayName', lbl);
  hold on
end
hold off
set(gca,'tickdir', 'out', 'fontsize',12)
set(gca,'fontsize',11)
xlabel('$t$', 'interpreter', 'latex', 'fontsize',12)
ylabel('$U(t)$', 'interpreter','latex', 'fontsize', 12, ...
  'rotation', 0)
legend('interpreter', 'latex','fontsize',12,'location','northeast',...
  'numcolumns',2)
ylim([-325 350])
%% state plot (Fig. 3)
f2 = figure(2);
f2.Position = 40*[5 5 15 13];
for k=1:4
  subplot(2,2,k)
  pn = nn(k); % state component to draw
  surf(TT, x1, cu{k}((pn-1)*xg+1:pn*xg,:))
  shading interp
  colormap winter
  view([-155 30])
  xlabel('$t$','interpreter', 'latex','fontsize', 12)
  ylabel('$x$','interpreter', 'latex','fontsize', 12)
  set(gca,'ytick',[0 1], 'fontsize',11)
  zl = ['$u^{', num2str(nn(k)), '}(t,x)$'];
  zlabel(zl,'interpreter', 'latex','fontsize', 12)
  if floor(k/2) == k/2
    set(gca,'position',get(gca,'position')+[.05 0 0 0])
  end
end