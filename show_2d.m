close all
%%

%vdy=zeros(nx0,ny0);
iz=2;%2d

%iz=6;%balloon
for nt=1:nts
dx=alx/(nx-1);
dy=aly/(ny-1);
xx=0:dx:alx;
yy=0:dy:aly;    
    load(['dat',sprintf('%4.4d',nt)])
 figure
    subplot(322)
    phi(:,:,iz)=phi(:,:,iz)+3;
    pcolor(yy,xx,squeeze(phi(:,:,iz))); 
    colorbar; shading interp;
    axis([0 100 10 65]);
  hold on
    ysep=40;
    lx=[0,100];ly=[ysep,ysep];
    plot(lx,ly,'k');
  hold off
    title('$$\phi$$','Interpreter','latex')
    ylabel('x')
    %xlabel('y')
    drawnow
    
    subplot(323)
    a=sp0(phi);
    pcolor(yy,xx,squeeze(a(:,:,2)))
    colorbar;shading interp;
    title('$$\widetilde{\phi}$$','Interpreter','latex')
    ylabel('x')
  hold on
    ysep=40;
    lx=[0,100];ly=[ysep,ysep];
    plot(lx,ly,'k');
  hold off
    axis([0 100 10 65]);
    %xlabel('y')
    drawnow
        
    subplot(324)
    pcolor(yy,xx,squeeze(wi(:,:,iz))); colorbar; shading interp;
    title('w')
  hold on
    ysep=40;
    lx=[0,100];ly=[ysep,ysep];
    plot(lx,ly,'k');
  hold off
    axis([0 100 10 65]);
    ylabel('x')
    %xlabel('y')
    drawnow

    % subplot(424)
    % a=sp0(pei);
    % pcolor(yy,xx,squeeze(a(:,:,iz)));colorbar; shading interp;
    % title('p')
    % xlabel('y')
    % ylabel('x')
    % drawnow

    subplot(321)
    pcolor(yy,xx,squeeze(pei(:,:,iz))); 
    caxis([0.2 1]);
    colorbar; shading interp;
    title('p_e')
    hold on
    ysep=40;
    lx=[0,100];ly=[ysep,ysep];
    plot(lx,ly,'k');
  hold off
    axis([0 100 10 65]);
    ylabel('x')
    %xlabel('y')
    drawnow

    subplot(325) 
   % a=sp0(pii);
     
    %vex=-1*dy2*(phi(2:nx-1,3:ny,iz)-phi(2:nx-1,1:ny-2,iz)); pcolor(vex);

    pcolor(yy,xx,squeeze(vex(:,:,iz))); 
    colorbar; shading interp; 
    title('V_x')
    hold on
    ysep=40;
    lx=[0,100];ly=[ysep,ysep];
    plot(lx,ly,'k');
  hold off
    axis([0 100 10 65]);
    xlabel('y')
    ylabel('x')
    drawnow

%     subplot(427)
%    % pcolor(squeeze(phi2(:,:,2))); colorbar; shading interp; title('phi2')
%    dx=alx/(nx-3);
%    dy=aly/(ny-3);
% xx=0:dx:alx;
% yy=0:dy:aly;
%    vdy=dx2*(pei(3:nx,2:ny-1,iz)-pei(1:nx-2,2:ny-1,iz));pcolor(yy,xx,vdy); colorbar; shading interp; title('dpdx')
  
%     xlabel('y')
%     ylabel('x')
%     drawnow

    subplot(326) 
    dx=alx/(nx-12);
dy=aly/(ny-1);
xx=0:dx:alx;
yy=0:dy:aly;
    %a=phi;
    %vey=dx2*(a(3:nx,2:ny-1,iz)-a(1:nx-2,2:ny-1,iz));pcolor(vey);
    pcolor(yy,xx,squeeze(vey(2:nx-10,:,iz))); 
    colorbar; shading interp; 
    axis([0 100 10 65]);
    title('V_y')
    hold on
    ysep=40;
    lx=[0,100];ly=[ysep,ysep];
    plot(lx,ly,'k');
  hold off
    xlabel('y')
    ylabel('x')
    drawnow
print(gcf,'-depsc',sprintf('%4.4d',nt))
close
end

%%


%    figure;
figure;plot(squeeze(wi(:,20,2)),'-o');
figure;plot(squeeze(phi(:,20,2)),'-o');
%%
figure;plot(squeeze(phi2(:,20,2)),'-o')
print(gcf,'-dpng','aphi2') 
%%
figure;plot(squeeze(phi(:,20,2)),'-o')
print(gcf,'-dpng','aphi') 
figure;plot(squeeze(pei(:,20,2)),'-o')
print(gcf,'-dpng','ape')
figure;plot(squeeze(wi(:,20,2)),'-o')
print(gcf,'-dpng','aw')
%Phi0=pei+phi;
%figure;plot(squeeze(Phi0(:,20,2)),'-o')
%print(gcf,'-dpng','aPhi0')
%figure;plot(squeeze(laplace2_pi(:,120,2)),'-o');print(gcf,'-dpng','laplace2pi')  
figure;plot(squeeze(vex(:,20,2)),'-o');
print(gcf,'-dpng','avey') 
%figure;
close all
%%
iz=2;
nts=763;
st=1;
for nt=st:nts
    
    load(['dat',sprintf('%4.4d',nt)])
    
    ptilde=std(pei(:,:,iz),0,2);
   %plot(squeeze(ptilde),'-o');print(gcf,'-dpng',sprintf('ptilde%4.4d',nt)) 
   
   p0=mean(pei(:,:,iz),2);
  plot(squeeze(ptilde./p0),'-o');
 % xlabel('x/\rho_s');ylabel('$\tilde{p}/p_0$','interpreter','latex');
  %print(gcf,'-depsc',sprintf('p0tilde%4.4d',nt)) 
  print(gcf,'-dpng',sprintf('p0tilde%4.4d',nt)) 
    
%plot(squeeze(p0),'-o');print(gcf,'-dpng',sprintf('p%4.4d',nt)) 
 
%plot(squeeze(mean(phi(:,:,2),2)),'-o');print(gcf,'-dpng',sprintf('phi%4.4d',nt)) 
%plot(squeeze(mean(vey(:,:,2),2)),'-o');print(gcf,'-dpng',sprintf('vy%4.4d',nt)) 
%plot(squeeze(mean(wi(:,:,2),2)),'-o');print(gcf,'-dpng',sprintf('w%4.4d',nt)) 

end
close all
%%
nts=313;
iz=2;
time=zeros(nts,1);time2=time;
time3=time;
nx=nx0+2;
ny=ny0+2;
dx=alx/(nx-1);
dx2=1./(2.*dx);
vdy=zeros(nx0,ny0);
x=0;
for nt=1:nts
	if nt<14
      x=100*0.01+x;
      t(nt)=x;
    elseif nt<114
    	x=50*0.005+x;
    	t(nt)=x;
    else
    	x=80*0.005+x;
    	t(nt)=x;
    end
end
ixs=25;

for nt=1:nts 
    load(['dat',sprintf('%4.4d',nt)])
    
vy2=mean(vey(:,:,iz).^2,2);
vy2=sqrt(vy2);
%plot(vy2,'-o');print(gcf,'-dpng',sprintf('rms_uy%4.4d',nt)) 

time(nt)=mean(vy2(ixs:nx2));

%
vdy=dx2*(pei(3:nx,2:ny-1,iz)-pei(1:nx-2,2:ny-1,iz));

%viy=vdy+squeeze(vey(2:nx0+1,2:ny0+1,iz));
%viy2=sqrt(mean(viy.^2,2));
%time3(nt)=mean(viy2(ixs:nx2));


vdy2=sqrt(mean(vdy.^2,2));
%plot(vdy,'-o');print(gcf,'-dpng',sprintf('rms_vdy%4.4d',nt)) 
time2(nt)=mean(vdy2(ixs:nx0));
end
%plot(time)
%hold 
plot(t,time2,'-r','Linewidth',1)

%plot(time3,'-m','Linewidth',1)
xlabel('t(R/c_s)');
legend('vdy');
%legend('vey','vdy','viy','Location','best');
print(gcf,'-dpng',sprintf('rms_uy')) 
%hold off
close all
%%
iz=2;
nts=313;
time=zeros(nts,1);time2=time;

for nt=1:nts   
    load(['dat',sprintf('%4.4d',nt)])
    %uxuy=vex(:,:,iz).*vey(:,:,iz);uxuytime=mean(uxuy,2).*mean(wi(:,:,iz),2);
    %plot(squeeze(uxuytime),'-o');print(gcf,'-dpng',sprintf('uxuy%4.4d',nt)) 

%time(nt)=mean(uxuytime);

vxsave=pei(90:102,:,iz).*vex(90:102,:,iz);
%vxsave=pei(:,:,iz).*vex(:,:,iz);
%plot(squeeze(mean(vxsave,2)),'-o');print(gcf,'-dpng',sprintf('flux%4.4d',nt)) 

time2(nt)=mean(mean(vxsave,2),1);

end
%time=100*time;
%plot(time)
%hold on
plot(time2,'-r','Linewidth',1)
xlabel('t'); legend('F_p');
%hold off
print(gcf,'-dpng',sprintf('Fp')) 

close all
    
%%
ny=ny0+2;
nts=413;
time=zeros(nts,1);
time2=time;
st=1;

for nt=st:nts
   load(['dat',sprintf('%4.4d',nt)])
   vy0=mean(vey(26:166,:,2),2);
   time(nt)=mean(vy0.^2);
   
   vy=vey(26:166,:,2)-repmat(vy0,[1,ny]);
   vperp=vy.^2+squeeze(vex(26:166,:,2).^2);
  time2(nt)=squeeze(mean(mean(vperp,2),1));
  
 
  
end
plot(time)
hold 
plot(time2,'-r','Linewidth',1)
xlabel('t'); legend('v0^2','v^2','Location','best');
print(gcf,'-dpng',sprintf('KE')) 
hold off
close all
%%
iz=2;
time=zeros(nts,1);
time2=time;
st=1;
for nt=st:nts
   load(['dat',sprintf('%4.4d',nt)])
   
edge=mean(pei(2:nx2,:,iz),2);
time(nt)=mean(edge);

sol=mean(pei(1+nx2:nx0+1,:,iz),2);
time2(nt)=mean(sol);
end
plot(time)
hold 
plot(time2,'-r','Linewidth',1)
xlabel('t'); legend('edge','sol','Location','best');
print(gcf,'-dpng',sprintf('psol')) 
hold off
close all
%%
%the confinement time
nt=20
nt=1
 load(['dat',sprintf('%4.4d',nt)])
nz=2;%2d
nxs=2*nx0*xw/alx;
nxs=round(nxs)
nx2
pro_pe=mean(pei(nxs:nx2,:,nz),2);
peAVG=sum(pro_pe)
SpAVG=sum(s(2:nx)/tau)
tau_cnfn=peAVG/SpAVG


