global ny nx

dx=alx/(nx-1);
for i=1:nx 
   x(i)=dx*(i-1); 
end 


%% The average profile 
st=240;
nts=260;
for nt=st:nts
   load(['dat',sprintf('%4.4d',nt)])
   phisave(:,:,nt)=phi(:,:,2);
end

   time3=mean(mean(phisave(:,:,nt),2),3); 

figure; 
   set(gca,'FontSize',14);
   plot(x,time3,'-','Linewidth',1)
   xlabel('x/\rho_s'); 
   legend('\phi');
   print(gcf,'-dpng','profphi'); 

%%
st=35;
nts=100;
for nt=st:nts
   load(['dat',sprintf('%4.4d',nt)])
   pesave(:,:,nt)=pei(:,:,2);
end
   time2=mean(mean(pesave(:,:,nt),2),3); 
figure; 
   set(gca,'FontSize',14);
   plot(x,time2,'r--','Linewidth',1)    
   xlabel('x/\rho_s'); 
   legend('P');
   print(gcf,'-dpng','profden_P'); 


%% Flux
st=1;
nts=60;
for nt=st:nts
   load(['dat',sprintf('%4.4d',nt)])
   vxsave2(:,:,nt)=pei(:,:,2).*vex(:,:,2);
end
time2=mean(mean(vxsave2(:,:,nt),2),3); 

figure; 
   set(gca,'FontSize',14);
   plot(x,time2,'r--','Linewidth',1)
   xlabel('x/\rho_s'); 
   ylabel('\Gamma(p_0^{.}c_s)');
   legend('p^{.}v_x')
   print(gcf,'-dpng','profflux')


%% Vy and W
for nt=st:nts
   load(['dat',sprintf('%4.4d',nt)])
   vysave(:,:,nt)=vey(:,:,2);
   wsave(:,:,nt)=wi(:,:,2); 
end
time=mean(mean(vysave(:,:,nt),2),3); 
time2=mean(mean(wsave(:,:,nt),2),3); 
figure; 
   set(gca,'FontSize',14);
   plot(x,time,'r--','Linewidth',1)
hold; 
   plot(x,time2,'m--','Linewidth',1) 
   xlabel('x/\rho_s'); 
   legend('v_y','w');
hold off;
   print(gcf,'-dpng','profw')


