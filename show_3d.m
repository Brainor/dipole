close all

figure;plot( cx,'-o');print(gcf,'-dpng',['cx','.png']); close
figure;plot( cy,'-o');print(gcf,'-dpng',['cy','.png']);close
%%

%pause
for nt=1:nts
figure; %nt
 load(['dat',sprintf('%4.4d',nt)])
     
 %phipi=pii+phi;ve=vii-jz./deni;
 %iz=floor(nz/2); 
%iz=nz/4+1;iz=(nz);
iz=(nz)/2+1;iz=10;%outside
iz=6;%outside
%iz=2;%2d
%iz=10;


   subplot(3,3,1)
    pcolor(pei(:,:,iz) ); colorbar;shading interp; 
title('p'); xlabel('y'); ylabel('x')
%axis image; 
drawnow
        subplot(3,3,2)
    pcolor(phi(:,:,iz)-0 ); colorbar;shading interp;
%axis image; title('\Phi')
title('potential'); xlabel('y'); ylabel('x')
drawnow
        subplot(3,3,3)
    pcolor(wi(:,:,iz) ); colorbar;shading interp;
title('vorticity'); xlabel('y'); ylabel('x')
%axis image; 
drawnow 
       subplot(3,3,4)
 a=pei; petilde=sp0(a);  a=squeeze(petilde(:,:,iz));  pcolor(a);
  %pcolor(pei(:,:,iz)); 
  colorbar; shading interp
title('p_e'); xlabel('y'); ylabel('x')
%axis image; 
drawnow
        subplot(3,3,5)
        
        %a=petilde+phi;
        %a=phi2;a=squeeze(a(:,:,iz));  pcolor(a);

 %       a=vii-jz;  a=squeeze(a(:,:,iz));  pcolor(a);
    %    a=vii; a=sp0(a);  a=squeeze(a(:,:,iz));  pcolor(a);
    pcolor(vii(:,:,iz));
    colorbar;shading interp
%title('u_{e}'); 
title('u_{i}');
%title('phi2');

xlabel('y'); ylabel('x')
%axis image; 
drawnow
        subplot(3,3,6)
    pcolor(jz(:,:,iz)); 
    colorbar;shading interp
 title('J_z'); 
  xlabel('y'); ylabel('x')        
%axis image;
drawnow 
        subplot(3,3,7)
 a=vex(:,:,iz); a=squeeze(a); pcolor(a); colorbar; shading interp
title('v_x'); xlabel('y'); ylabel('x')
       

        subplot(3,3,8)
                
a=vey(2:nx0-10,:,iz); a=squeeze(a); pcolor(a); colorbar; shading interp
title('v_y');

xlabel('y'); ylabel('x')
        subplot(3,3,9)
a=phi; a=sp0(a);a=squeeze(a(:,:,iz)); pcolor(a); colorbar; shading interp
title('phi'); xlabel('y'); ylabel('x')

print(gcf,'-dpng',sprintf('%2.2dz%4.4d',iz,nt)) 
close 
end 

%%

%%
i2=ny/2;
for nt=1:nts
figure; 
 
load(['dat',sprintf('%4.4d',nt)])

    subplot(3,3,1)
a=pei(:,i2,:); a=squeeze(a); pcolor(a); colorbar; shading interp
title('p'); xlabel('z'); ylabel('x')
        subplot(3,3,2)
a=phi-0;a=squeeze(a(:,i2,:)); pcolor(a); colorbar; shading interp
title('potential'); xlabel('z'); ylabel('x')
        subplot(3,3,3)
a=wi(:,i2,:); a=squeeze(a); pcolor(a); colorbar; shading interp
title('vorticity'); xlabel('z'); ylabel('x')
    subplot(3,3,4)
%a=pei(:,i2,:); 
a=pei; a=sp0(a);  a=squeeze(a(:,i2,:)); 
pcolor(a); colorbar; shading interp 
title('p_e'); xlabel('z'); ylabel('x')
        subplot(3,3,5)
a=vii(:,i2,:); a=squeeze(a); pcolor(a); colorbar; shading interp
title('u_i'); xlabel('z'); ylabel('x')
        subplot(3,3,6)
a=jz(:,i2,2:nz-1); a=squeeze(a); pcolor(a); colorbar;shading interp; title('J_z') ; xlabel('z'); ylabel('x')
        subplot(3,3,7)
        a=vex(:,i2,:); a=squeeze(a); pcolor(a); colorbar; shading interp
title('v_x'); xlabel('z'); ylabel('x')

        subplot(3,3,8)
          %a=vii-jz;  
         % a=vii;a=squeeze(a(1:nx2-5,i2,:));  pcolor(a);
%a=vey(:,i2,:); a=squeeze(a);  
a=wi; a=sp0(a);  a=squeeze(a(:,i2,:));
pcolor(a);
colorbar; shading interp
title('w'); xlabel('z'); ylabel('x')
        subplot(3,3,9)
a=phi; a=sp0(a);a=squeeze(a(:,i2,:)); pcolor(a); colorbar; shading interp
title('phi'); xlabel('z'); ylabel('x')

print(gcf,'-dpng',sprintf('%2.2dy%4.4d',i2,nt))
close
end


%%

%pause
for nt=1:nts
figure; %nt
 load(['dat',sprintf('%4.4d',nt)])
     
 %phipi=pii+phi;ve=vii-jz./deni;
 %iz=floor(nz/2); 
%iz=nz/4+1;iz=(nz);
iz=(nz)/2+1;iz=10;%outside
iz=6;%outside
%iz=2;%2d
%iz=10;


   subplot(3,3,1)
    pcolor(pei(:,:,iz) ); colorbar;shading interp; 
title('p'); xlabel('y'); ylabel('x')
%axis image; 
drawnow
        subplot(3,3,2)
    pcolor(phi(:,:,iz)-0 ); colorbar;shading interp;
%axis image; title('\Phi')
title('potential'); xlabel('y'); ylabel('x')
drawnow
        subplot(3,3,3)
    pcolor(wi(:,:,iz) ); colorbar;shading interp;
title('vorticity'); xlabel('y'); ylabel('x')
%axis image; 
drawnow 
       subplot(3,3,4)
 a=pei; a=sp0(a);  a=squeeze(a(:,:,iz));  pcolor(a);
  %pcolor(pei(:,:,iz)); 
  colorbar; shading interp
title('p_e'); xlabel('y'); ylabel('x')
%axis image; 
drawnow
        subplot(3,3,5)
        a=vii; a=sp0(a);  a=squeeze(a(:,:,iz));  pcolor(a);
    %pcolor(vii(:,:,iz));
    colorbar;shading interp
title('u_{i}'); xlabel('y'); ylabel('x')
%axis image; 
drawnow
        subplot(3,3,6)
    pcolor(jz(:,:,iz)); colorbar;shading interp
 title('J_z'); xlabel('y'); ylabel('x')        
%axis image;
drawnow 
        subplot(3,3,7)
 a=vex(:,:,iz); a=squeeze(a); pcolor(a); colorbar; shading interp
title('v_x'); xlabel('y'); ylabel('x')
       

        subplot(3,3,8)
                a=vii-jz;  a=squeeze(a(:,:,iz));  pcolor(a);
    %pcolor(vii(:,:,iz));
    colorbar;shading interp
title('u_{e}'); xlabel('y'); ylabel('x')
%a=vey(:,:,iz); a=squeeze(a); pcolor(a); colorbar; shading interp
%title('v_y');

xlabel('y'); ylabel('x')
        subplot(3,3,9)
a=phi; a=sp0(a);a=squeeze(a(:,:,iz)); pcolor(a); colorbar; shading interp
title('phi'); xlabel('y'); ylabel('x')

print(gcf,'-dpng',sprintf('%2.2dz%4.4d',iz,nt)) 
close 
end 

%%

%pause
for nt=1:nts
figure; %nt
load(['dat',sprintf('%4.4d',nt)])

 %phipi=pii+phi;ve=vii-jz./deni;
 %iz=floor(nz/2); 
%iz=nz/4+1;iz=(nz);
iz=2;%inside
iz=14;%inside
%iz=(nz)/2+1;


   subplot(4,2,1)
    pcolor(deni(:,:,iz) ); colorbar;shading interp; 
title('n'); xlabel('y'); ylabel('x')
%axis image; 
drawnow
        subplot(4,2,2)
    pcolor(phi(:,:,iz)-0 ); colorbar;shading interp;
%axis image; title('\Phi')
title('potential'); xlabel('y'); ylabel('x')
drawnow
        subplot(4,2,3)
    pcolor(wi(:,:,iz) ); colorbar;shading interp;
title('vorticity'); xlabel('y'); ylabel('x')
%axis image; 
drawnow 
    subplot(4,2,4)
 a=pei; a=sp0(a);  a=squeeze(a(:,:,iz));  pcolor(a);
  %pcolor(pei(:,:,iz)); 
  colorbar; shading interp
title('p_e'); xlabel('y'); ylabel('x')
%axis image; 
drawnow
        subplot(4,2,5)
    pcolor(vii(:,:,iz)); colorbar;shading interp
title('u_{i}'); xlabel('y'); ylabel('x')
%axis image; 
drawnow
        subplot(4,2,6)
    pcolor(jz(:,:,iz)); colorbar;shading interp
 title('J_z'); xlabel('y'); ylabel('x')
%axis image;
drawnow 
        subplot(4,2,7)
a=vex(:,:,iz); a=squeeze(a); pcolor(a); colorbar; shading interp
title('v_x'); xlabel('y'); ylabel('x')
        subplot(4,2,8)
a=vey(:,:,iz); a=squeeze(a); pcolor(a); colorbar; shading interp
title('v_y'); xlabel('y'); ylabel('x')

print(gcf,'-dpng',sprintf('%2.2dz%4.4d',iz,nt)) 
close 
end 
%%
%i2=25;
i2=10;
%plot y-z

for nt=1:nts
 figure;%nt
 load(['dat',sprintf('%4.4d',nt)])

%phipi=pii+phi;ve=vii-jz./deni;
    subplot(3,3,1)
a=deni(i2,:,:); a=squeeze(a); pcolor(a); colorbar; shading interp
title('n'); xlabel('z'); ylabel('y')
        subplot(3,3,2)
a=phi-0; a=squeeze(a(i2,:,:)); pcolor(a); colorbar; shading interp
title('potential'); xlabel('z'); ylabel('y')
        subplot(3,3,3)
a=wi(i2,:,:); a=squeeze(a); pcolor(a); colorbar; shading interp
title('vorticity'); xlabel('z'); ylabel('y')
    subplot(3,3,4)
a=pei;%a=sp0(a);  
a=squeeze(a(i2,:,:)); pcolor(a); colorbar; shading interp 
title('p_e'); xlabel('z'); ylabel('y')
        subplot(3,3,5)
a=vii(i2,:,:); a=squeeze(a); pcolor(a); colorbar; shading interp
title('u_{i}'); xlabel('z'); ylabel('y')
        subplot(3,3,6)
a=jz(i2,:,:); a=squeeze(a); pcolor(a); colorbar;shading interp; 
title('J_z') ; xlabel('z'); ylabel('y')
        subplot(3,3,7)
a=pii; a=squeeze(a(i2,:,:)); pcolor(a); colorbar; shading interp
title('P_i'); xlabel('z'); ylabel('y')
        subplot(3,3,8)
a=vex(i2,:,:); a=squeeze(a); pcolor(a); colorbar; shading interp
title('v_x'); xlabel('z'); ylabel('y')
        subplot(3,3,9)
a=vey(i2,:,:); a=squeeze(a); pcolor(a); colorbar; shading interp
title('v_y'); xlabel('z'); ylabel('y')

print(gcf,'-dpng',sprintf('%2.2dx%4.4d',i2,nt))
close
end

%%
%%
for nt=1:nts
figure; %nt
 
load(['dat',sprintf('%4.4d',nt)])

 %phipi=pii+phi;ve=vii-jz./deni;
 %iz=floor(nz/2); 
%a=phi-0;
a=sp0(phi);
%a=wi;
%a=pei;
%a=den0;a=deni;
iz=1;

   subplot(2,3,1)
    pcolor(a(:,:,iz) ); colorbar;shading interp; 
title(sprintf('z%2.2d',iz)); xlabel('y'); ylabel('x')
%axis image; 
drawnow
iz=2;

   subplot(2,3,2)
    pcolor(a(:,:,iz) ); colorbar;shading interp; 
title(sprintf('z%2.2d',iz)); xlabel('y'); ylabel('x')
%axis image; 
drawnow

iz=6;
        subplot(2,3,3)
    pcolor(a(:,:,iz) ); colorbar;shading interp;
%axis image; title('\Phi')
title(sprintf('z%2.2d',iz)); xlabel('y'); ylabel('x')
drawnow
%iz=14;
iz=10;
        subplot(2,3,4)
    pcolor(a(:,:,iz) ); colorbar;shading interp;
title(sprintf('z%2.2d',iz)); xlabel('y'); ylabel('x')
%axis image; 
drawnow 
iz=nz-1;
    subplot(2,3,5)
   pcolor(a(:,:,iz));
  %pcolor(pei(:,:,iz)); 
  colorbar; shading interp
%title('');
title(sprintf('z%2.2d',iz));xlabel('y'); ylabel('x')
%axis image; 
drawnow
iz=nz;
    subplot(2,3,6)
  a=squeeze(a(:,:,iz));  pcolor(a);
  %pcolor(pei(:,:,iz)); 
  colorbar; shading interp
%title('');
title(sprintf('z%2.2d',iz));xlabel('y'); ylabel('x')
%axis image; 
drawnow
  

print(gcf,'-dpng',sprintf('phi_tilde%4.4d',nt)) 
close 
end 

%%

%figure;plot(squeeze(pei(10,:,10)),'-o')
%print(gcf,'-dpng',['deny','.png'])
%figure;plot(squeeze(phi(10,:,10)),'-o')
%print(gcf,'-dpng',['phiy','.png'])
%%
a=sp0(pei);
figure;plot(squeeze(a(:,10,2)),'-o');
print(gcf,'-dpng',['z2pertx','.png'])

figure;plot(squeeze(a(10,10,:)),'-o');
print(gcf,'-dpng',['pertz','.png'])

figure;plot(squeeze(pei(10,10,:)),'-o');
print(gcf,'-dpng',['pez','.png'])
%%
figure;plot(squeeze(vii(10,10,:)),'-o');print(gcf,'-dpng',['vz','.png'])
figure;plot(squeeze(tei(10,10,:)),'-o');
print(gcf,'-dpng',['tz','.png'])
%%
figure;plot(squeeze(wi(10,10,:)),'-o');
print(gcf,'-dpng',['wz','.png'])
%%
figure;plot(squeeze(jz(10,10,:)),'-o');
print(gcf,'-dpng',['jz','.png'])
%%
figure;plot(squeeze(phi(10,10,:)),'-o')
print(gcf,'-dpng',['phiz','.png'])
%%
figure;plot(squeeze(phi(:,10,2)),'-o');
print(gcf,'-dpng',['z2phix','.png'])
%%
figure;plot(squeeze(wi(:,10,1)),'-o');
print(gcf,'-dpng',['z1wx','.png'])
%%
figure;plot(squeeze(pei(:,10,2)),'-o');
print(gcf,'-dpng',['z2pex','.png'])
%%
figure;plot(squeeze(vii(:,10,2)),'-o');
print(gcf,'-dpng',['z2ux','.png'])
%%
figure;plot(squeeze(jz(:,10,2)),'-o');
print(gcf,'-dpng',['z2jx','.png'])
%%
%plot open
figure;plot(squeeze(pei(40,10,:)),'-o');
print(gcf,'-dpng',['x40denz','.png'])
%%
figure;plot(squeeze(vii(40,10,:)),'-o');
print(gcf,'-dpng',['x40vz','.png'])
figure;plot(squeeze(tei(40,10,:)),'-o');
print(gcf,'-dpng',['x40t','.png'])
figure;plot(squeeze(wi(40,10,:)),'-o');
print(gcf,'-dpng',['x40w','.png'])
%%
figure;plot(squeeze(jz(40,10,2:nz-1)),'-o');
print(gcf,'-dpng',['x40jz','.png'])

figure;plot(squeeze(phi(40,10,:)),'-o')
print(gcf,'-dpng',['x40phiz','.png'])


close all

%figure;plot(squeeze(pei(10,10,:)),'-o')
%figure;plot(squeeze(pii(10,10,:)),'-o')
%figure;plot(squeeze(vii(10,10,:)),'-o')
%figure;plot(squeeze(jz(10,10,:)),'-o')


