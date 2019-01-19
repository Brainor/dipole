% use this in steady state, analyze < u_x u_y >
% maybe Reynolds stress
% choose specific z, mean in y and time
% changing in x
close all
%clear all
nx=nx0+2;
ny=ny0+2;
nx2=102;
dx=alx/(nx-1); % fixed B.C in inner part @x
dy=aly/(ny0);

% na, nb denotes the time range
na=1;
nb=100;
iz=2;

% pettime means tilde of pe in time
uxuytime=zeros(nx,nb-na+1);
tn=0;
for nt=na:nb
    nt
    tn=tn+1;
    load(['dat',sprintf('%4.4d',nt)])
    uxuy=vex(:,:,iz).*vey(:,:,iz);
    uxuytime(:,tn)=mean(uxuy,2).*mean(wi(:,:,2),2);

end

% time average

xx=0:dx:alx;

uxuym=mean(uxuytime,2);
ym=max(uxuym);

figure;
% set(gca,'fontsize',14);
plot(xx,uxuym)


hold on
% plot(xx,p0,'b:');
xsep=nx2/nx*alx;
lx=[xsep,xsep];ly=[-5*10^-5,1.1*ym];
plot(lx,ly,'k--');
text(xsep+1,1.1*ym,'LCFS','fontsize',14);
get(gca);
xlabel('x/\rho_s');
ylabel('< u_x u_y>w');

drawnow;
print(gcf,'-dpng',sprintf('uxuywt%2.2dz%2.2d-%2.2d',iz,na,nb));
close

