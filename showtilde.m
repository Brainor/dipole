% use this in steady state
% choose specific z, mean in y and time
% changing in x
close all
%clear all
nx=nx0+2;
ny=ny0+2;
dx=alx/(nx-1); % fixed B.C in inner part @x
dy=aly/(ny0);

% na, nb denotes the time range
na=10;
nb=30;
iz=2;
% pettime means tilde of pe in time
pettime=zeros(nx,nb-na+1);
denttime=pettime;
pers=pettime;
denrs=pettime;
p0time=pettime;
for nt=na:nb
    nt
    
    load(['dat',sprintf('%4.4d',nt)])
    
    p0time(:,nt)=mean(pei(:,:,iz),2);
    pettime(:,nt)=std(pei(:,:,iz),1,2);
    pers(:,nt)=pettime(:,nt)./p0time(:,nt);
    denttime(:,nt)=std(deni(:,:,iz),1,2);
    denrs(:,nt)=denttime(:,nt)./mean(deni(:,:,iz),2);
end

% time average
pet=mean(pettime,2);
% dent=mean(denttime,2);
p0=mean(p0time,2);
pert=mean(pers,2);
% denrt=mean(denrs,2);

%     pet=100*pet;

xx=0:dx:alx;
figure;
ym=max(pet);
% set(gca,'fontsize',14);
[AX,H1,H2]=plotyy(xx,pet,xx,p0,'plot');
legend('pet','p_0');
set(get(AX(1),'Ylabel'),'String','pet','fontsize',14) 
set(get(AX(2),'Ylabel'),'String','p_0')%,'fontsize',14)

set(H1,'LineStyle','-');
set(H2,'Marker','+');
hold on
% plot(xx,p0,'b:');
xsep=1/2*alx;
lx=[xsep,xsep];ly=[0,1.1*ym];
plot(lx,ly,'k--');
text(xsep+1,1.1*ym,'LCFS','fontsize',14);
get(gca);
xlabel('x/\rho_s');
ylabel('$\tilde{p}_e$','interpreter','latex');

drawnow;
print(gcf,'-dpng',sprintf('testpet%2.2dz%2.2d-%2.2d',iz,na,nb));
close

% ym=max(pet);
% set(gca,'fontsize',14);
% plotyy(xx,pet,'r-',xx,p0,'b+');
% 
% hold on
% % plot(xx,p0,'b:');
% xsep=nx2/nx*alx;
% lx=[xsep,xsep];ly=[0,1.1*ym];
% plot(lx,ly,'k--');
% text(xsep+1,1.1*ym,'LCFS','fontsize',14);
% get(gca);
% xlabel('x/\rho_s');
% ylabel('$\tilde{p}_e$','interpreter','latex');
% legend('pet','p0','LCFS');
% drawnow;
% print(gcf,'-depsc',sprintf('pet%2.2dz%2.2d-%2.2d',iz,na,nb));
% close

figure;
ym=max(pert);
set(gca,'fontsize',14);
plot(xx,pert,'r-');
hold on
xsep=1/2*alx;
lx=[xsep,xsep];ly=[0,1.2*ym];
plot(lx,ly,'k--');
text(xsep+1,1.2*ym,'LCFS','fontsize',14);
%axis([0 25 0 0.18])
get(gca);
xlabel('x/\rho_s');
ylabel('$\tilde{p}_e/p_0$','interpreter','latex');
drawnow;
print(gcf,'-dpng',sprintf('pert%2.2dz%2.2d-%2.2d',iz,na,nb));
close

% gpe=firstlog(pet,dx);

% plot(gpe(1:30))
% set(gca,'fontsize',14)
% xlabel('t')
% ylabel('\gamma')
% title('$\gamma of \tilde{p}_e -t$','interpreter','latex');
% print(gcf,'-deps','gamma');
