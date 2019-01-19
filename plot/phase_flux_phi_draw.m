close all
% alx=75;nx0=192;
% aly=100;ny0=256;
% nx=nx0+2;
% ny=ny0+2;

% ix=77;
% iz=2;

% nk=2;
% st=190;
% nts=1190;
% x=0;
% for nt=1:nts
%   if nt<161
%       x=80*0.005+x;
%     else
%       x=100*0.002+x;
%     end
%     t(nt)=x;
% end
% petilde=zeros(ny0,1);
% vx=petilde;
% pet=zeros(nts-st+1,1);
% vxt=pet;
% pek=zeros(ny0,1);
% vxk=pek;
% phase=zeros(nts-st+1,1);
% cos_phase=phase;
% abs_flux=phase;
% vy=phase;

% for nt=st:nts
% load(['dat',sprintf('%4.4d',nt)])

% petilde(:)=pei(ix,2:ny0+1,iz)-mean(pei(ix,2:ny0+1,iz));
% vx(:)=vex(ix,2:ny0+1,iz);
% % phitilde(:)=phi(ix,2:ny0+1,iz)-mean(phi(ix,2:ny0+1,iz));

% pet(nt-st+1)=rms(petilde);
% vxt(nt-st+1)=rms(vx);

% pek(:)=ifft(petilde(:));
% vxk(:)=ifft(vx(:));

% phase(nt-st+1)=angle(pek(nk)/vxk(nk));
% cos_phase(nt-st+1)=cos(phase(nt-st+1));
% abs_flux(nt-st+1)=abs(pek(nk))*abs(vxk(nk));
% vy(nt-st+1)=mean(mean(vey(60:90,2:ny0+1,iz),1),2);

% end
nk=3;
load(['phase_flux_phi_nk3'])
figure;
height=0.19;
width=0.8;
left=0.12;
bottom=0.1;
ps1=[left,bottom,width,height];
ps2=[left,bottom+height+0.03,width,height];
ps3=[left,bottom+2*height+0.06,width,height];
ps4=[left,bottom+3*height+0.09,width,height];

subplot('position',ps4) 
set(gca,'FontSize',14);
plot(t(st:nts),phase/pi,'b-','Linewidth',1);
hold on;
xsep=110;
lx=[xsep,xsep];ly=[-1,1];
plot(lx,ly,'k--');
xsep=170;
lx=[xsep,xsep];ly=[-1,1];
plot(lx,ly,'k--');
ysep=0.5;
ly=[ysep,ysep];lx=[70,270];
plot(lx,ly,'r--');
ysep=-0.5;
ly=[ysep,ysep];lx=[70,270];
plot(lx,ly,'r--');
ylabel('\theta_k/\pi');
axis([70 270 -1 1]);
a=findall(gcf,'type','axes');
set(a,'XTickLabel',[]);
% xlabel('t (a/c_s)');

subplot('position',ps3) 
set(gca,'FontSize',14);
plot(t(st:nts),sin_phase,'b-','Linewidth',1);
hold on;
xsep=110;
lx=[xsep,xsep];ly=[-1,1];
plot(lx,ly,'k--');
xsep=170;
lx=[xsep,xsep];ly=[-1,1];
plot(lx,ly,'k--');
ysep=0;
ly=[ysep,ysep];lx=[70,270];
plot(lx,ly,'r--');
ylabel('sin\theta_k');
axis([70 270 -1 1]);
a=findall(gcf,'type','axes');
set(a,'XTickLabel',[]);

subplot('position',ps2) 
set(gca,'FontSize',14);
plot(t(st:nts),abs_flux,'b-','Linewidth',1);
hold on;
xsep=110;
lx=[xsep,xsep];ly=[0,1*10^-2];
plot(lx,ly,'k--');
xsep=170;
lx=[xsep,xsep];ly=[0,1*10^-2];
plot(lx,ly,'k--');
ylabel('k|p_k||\phi_k|');
axis([70 270 0 1*10^-2]);
a=findall(gcf,'type','axes');
set(a,'XTickLabel',[]);

subplot('position',ps1) 
set(gca,'FontSize',14);
plot(t(st:nts),vy,'b-','Linewidth',1);
hold on;
xsep=110;
lx=[xsep,xsep];ly=[0,0.4];
plot(lx,ly,'k--');
xsep=170;
lx=[xsep,xsep];ly=[0,0.4];
plot(lx,ly,'k--');
ylabel('v_y');
axis([70 270 0 0.4]);
xlabel('t (a/c_s)');

print(gcf,'-depsc',['phase_flux_phi_nk=',num2str(nk),'.eps']);
print(gcf,'-dpng',['phase_flux_phi_nk=',num2str(nk),'.png']);
% fid=['phase_flux_nk2'];
% save(fid,'t','phase','cos_phase','abs_flux','vy')
