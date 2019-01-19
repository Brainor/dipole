figure;
height=0.43;
width=0.89;
left=0.06;
bottom=0.08;
ps1=[left,bottom,width,height];
ps2=[left,bottom+height,width,height];

subplot('position',ps2)
plot(Time,K,'r-');
hold on
a=findall(gcf,'type','axes');
set(a,'XTickLabel',[]);
plot(Time,U,'b--');
% xlabel('$t(\rho_s/c_s)$','Interpreter','latex')
legend('K','U');
% title('kinetic energy and mean motions');
% text(9,36,'kinetic energy and mean motions');
drawnow

subplot('position',ps1)
plot(Time,peEdge,'r-');
hold on
plot(Time,peSOL,'b--');
xlabel('$t(\rho_s/c_s)$','Interpreter','latex')
legend('p_{edge}','p_{SOL}');
% title('pressure in edge and SOL');
% text(9,2.3,'pressure in edge and SOL');
drawnow

print(gcf,'-dpng','pES-KaU')