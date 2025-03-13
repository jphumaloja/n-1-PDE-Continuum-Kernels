%% Collect simulation data
nn = [3 5 10 20]; % selected values of n
% initialize
cu = cell(size(nn));
cku = cu;
eu = cu;
eku = cu;
% generic simulation parameters
T = 10; % simulation end time
xg = 256; % x grid
tmp = linspace(0,1,xg+1); % auxiliary grid
x0 = tmp(1:xg); x1 = tmp(2:xg+1); % grids for u and v
tg = 513; % time grid (for plotting the solution)
TT = linspace(0,T,tg); % time vector of solution
% simulate for fixed initial conditions (depending on n)
for k = 1:numel(nn)
  % simulate with continuum kernels
  [usol, Usol] = n1sim(nn(k),T,tg);
  cu{k} = usol;
  cku{k} = Usol;
  % and with exact kernel
  [usol, Usol] = n1sim(nn(k),T,tg,ksol(nn(k)));
  eu{k} = usol;
  eku{k} = Usol;
end
%% plot (Fig. 4)
f3 = figure(3);
f3.Position = 40*[5 5 15 12];
for k=1:numel(nn)
  subplot(numel(nn),1,k)
  plot(TT, [cku{k}; eku{k}],'-','linewidth',2)
  set(gca,'tickdir', 'out', 'fontsize',11)
  ylbl = ['$U_{',num2str(nn(k)),'}$'];
  ylabel(ylbl, 'interpreter','latex', 'fontsize', 12, ...
  'rotation', 0)
  legend({'continuum','exact'},'interpreter', 'latex',...
  'fontsize',12,'location','southeast','numcolumns',1)
  if k < numel(nn)
    set(gca,'xticklabel',{})
    set(gca,'position',get(gca,'position')+[0 0 .05 .1/numel(nn)])
  else
    set(gca,'position',get(gca,'position')+[0 0 .05 .1/numel(nn)])
    xlabel('$t$', 'interpreter', 'latex', 'fontsize',12)
  end
end
