
% 2d

global ny nx xmin xmax

dx=alx/(nx-1);
% for i=1:nx 
%    x(i)=dx*(i-1); 
% end 
x=xmin:dx:xmax;

%% The average profile of phi p
nts=200;
st=1;
for nt=st:nts
   load(['dat',sprintf('%4.4d',nt)])
   time3=mean(phi(:,:,2),2); 
  figure; 
   set(gca,'FontSize',14);
   plot(x,time3,'-','Linewidth',1)
hold;
   % a=alx*nx2/nx0;
   % y=-1.5:0.001:1.5;
   % plot(a,y,'-','Linewidth',1)
   % text(a+1,1.4,'LCFS','fontsize',14);
    xlabel('x/\rho_s'); 
    legend('\phi');
    % axis([0 75 -1.5 1.5]);
hold off;
   print(gcf,'-dpng',sprintf('profphi%4.4d',nt)) 
end
%%
st=1;
nts=50;
for nt=st:nts
   load(['dat',sprintf('%4.4d',nt)])
   time=mean(pei(:,:,2),2);
figure; 
   set(gca,'FontSize',14);
   plot(x,time,'-','Linewidth',1)
    xlabel('x/\rho_s'); 
    legend('p');
    axis([0 75 0.1 1.4]);
   print(gcf,'-dpng',sprintf('profpe%4.4d',nt)) 
end


%% Flux
st=90;
nts=210;
for nt=st:nts
   load(['dat',sprintf('%4.4d',nt)])
   %vxsave(:,:,nt)=deni(:,:,2).*vex(:,:,2);
   vxsave2(:,:,2)=pei(:,:,2).*vex(:,:,2);
%time=mean(mean(vxsave,2),3); 
time2=mean(vxsave2(:,:,2),2); 

figure; 
   set(gca,'FontSize',14);
   plot(x,time2,'-','Linewidth',1)

  
   xlabel('x/\rho_s'); 
   ylabel('\Gamma(p_0^{.}c_s)');
   legend('p^{.}v_x')
   %axis([0 75 -2*10^-3 4*10^-3]);

   % a=alx*nx2/nx0;
   % y=-1*10^-3:0.00001:10*10^-3;
   % plot(a,y,'-','Linewidth',1)
   % text(a+1,9*10^-3,'LCFS','fontsize',14);

   print(gcf,'-dpng',sprintf('profflux%4.4d',nt))
end


%% Vy
for nt=st:nts
   load(['dat',sprintf('%4.4d',nt)])

time=mean(vey(:,:,2),2); 

figure; 
   set(gca,'FontSize',14);
   plot(x,time,'r--','Linewidth',1)
hold;
   xlabel('x/\rho_s'); 
   legend('v_y(c_s)');
   % axis([0 65 -0.1 0.05]);

   % a=alx*nx2/nx0;
   % y=-0.1:0.00001:0.05;
   % plot(a,y,'-','Linewidth',1)
   % text(a+1,0.04,'LCFS','fontsize',14);

hold off;
   print(gcf,'-dpng',sprintf('profvy%4.4d',nt))
end


%% w

for nt=st:nts
   load(['dat',sprintf('%4.4d',nt)])
   % wsave(:,:,nt)=wi(:,:,2); 
time2=mean(wi(:,:,2),2); 
figure; 
   set(gca,'FontSize',14);
   plot(x,time2,'r--','Linewidth',1)
hold;    
   xlabel('x/\rho_s'); 
   legend('w(\Omega)');
   % axis([0 65 -0.02 0.02]);

   % a=alx*nx2/nx0;
   % y=-0.02:0.0001:0.02;
   % plot(a,y,'-','Linewidth',1)
   % text(a+1,0.015,'LCFS','fontsize',14);

hold off;
   print(gcf,'-dpng',sprintf('profw%4.4d',nt))
end