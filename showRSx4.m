
close all
%clear all
nx=nx0+2;
ny=ny0+2;
nx2=102;
dx=alx/(nx-1); % fixed B.C in inner part @x
dy=aly/(ny0);


na=124;
nb=134;
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
xx=0:dx:alx;
uxuym=mean(uxuytime,2);
ym=max(uxuym);
figure;
height=0.43;
width=0.85;
left=0.1;
bottom=0.1;
ps1=[left,bottom,width,height];
ps2=[left,bottom+height,width,height];
subplot('position',ps2) 
plot(xx,uxuym,'k')

hold on

na=174;
nb=184;
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
xx=0:dx:alx;
uxuym=mean(uxuytime,2);
ym=max(uxuym);
plot(xx,uxuym,'k--')

na=185;
nb=213;
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
xx=0:dx:alx;
uxuym=mean(uxuytime,2);
ym=max(uxuym);
plot(xx,uxuym,'k.')

xsep=40;
lx=[xsep,xsep];ly=[-8*10^-5,2*10^-4];
plot(lx,ly);
text(xsep+1,1.9*10^-4,'LCFS','fontsize',10);
get(gca);
a=findall(gcf,'type','axes');
set(a,'XTickLabel',[]);
ylabel('< u_x u_y>w');
legend('L','before transition','H');
axis([0 75 -8*10^-5 2*10^-4]);
hold off;




subplot('position',ps1)
na=124;
nb=134;
iz=2;
uxuytime=zeros(nx,nb-na+1);
tn=0;
for nt=na:nb
    nt
    tn=tn+1;
    load(['dat',sprintf('%4.4d',nt)])
    uxuy=vex(:,:,iz).*vey(:,:,iz);
    uxuytime(:,tn)=mean(uxuy,2);
end
xx=0:dx:alx;
uxuym=mean(uxuytime,2);
ym=max(uxuym);
plot(xx,uxuym,'k')

hold on
na=174;
nb=184;
iz=2;
uxuytime=zeros(nx,nb-na+1);
tn=0;
for nt=na:nb
    nt
    tn=tn+1;
    load(['dat',sprintf('%4.4d',nt)])
    uxuy=vex(:,:,iz).*vey(:,:,iz);
    uxuytime(:,tn)=mean(uxuy,2);
end
xx=0:dx:alx;
uxuym=mean(uxuytime,2);
ym=max(uxuym);
plot(xx,uxuym,'k--')

na=185;
nb=213;
iz=2;
uxuytime=zeros(nx,nb-na+1);
tn=0;
for nt=na:nb
    nt
    tn=tn+1;
    load(['dat',sprintf('%4.4d',nt)])
    uxuy=vex(:,:,iz).*vey(:,:,iz);
    uxuytime(:,tn)=mean(uxuy,2);

end
xx=0:dx:alx;
uxuym=mean(uxuytime,2);
ym=max(uxuym);
plot(xx,uxuym,'k.')

xsep=40;
lx=[xsep,xsep];ly=[-1*10^-3,7*10^-3];
plot(lx,ly);
text(xsep+1,6.5*10^-3,'LCFS','fontsize',10);
legend('L','before transition','H');
get(gca);
axis([0 75 -1*10^-3 7*10^-3]);
xlabel('x/\rho_s');
ylabel('< u_x u_y>');

drawnow;
print(gcf,'-depsc',sprintf('uxuyt%2.2dz',iz));