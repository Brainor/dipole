%»æÖÆÆµÆ×
close all
alx=75;nx0=192;
aly=100;ny0=256;
nx=nx0+2;
ny=ny0+2;
ntp=100;
tau=0.002;
ix=77; % x=30 rhos
iz=2;

st=191;         %% L
nts=390;
Nt=nts-st+1;
at=ntp*tau*Nt;
nt02=Nt/2;
      pi2=2.*pi;
      pi2t=pi2/at; 
     vf=zeros(Nt,1);
      for j=1:nt02+1
         vf(j)=pi2t*(j-1);
      end 
      for j=nt02+2:Nt
         jj=j-Nt;
         vf(j)=pi2t*(jj-1);
      end 
cr=zeros(Nt,1);
phif=zeros(Nt,1);
phit=phif;


for nt=st:nts
load(['dat',sprintf('%4.4d',nt)])
phit(nt-st+1)=mean(phi(ix,2:ny0+1,iz));
end
phit=phit-mean(phit); 
cr(:)=ifft(phit);
phif(:)=abs(cr(:));

figure; 
set(gca,'FontSize',14);
a(1)=phif(1);
a(2:Nt/2)=2*phif(2:Nt/2);

fid=['phif_L',sprintf('%4.4d',i)];
save(fid,'a')

plot(vf(1:Nt/2),a(1:Nt/2),'b-','Linewidth',1)
hold on;


st=401;         %% LCO
nts=500;
Nt=nts-st+1;
at=ntp*tau*Nt;
nt02=Nt/2;
pi2=2.*pi;
pi2t=pi2/at; 
vf=zeros(Nt,1);
      for j=1:nt02+1
         vf(j)=pi2t*(j-1);
      end 
      for j=nt02+2:Nt
         jj=j-Nt;
         vf(j)=pi2t*(jj-1);
      end 
cr=zeros(Nt,1);
phif=zeros(Nt,1);
phit=phif;

for nt=st:nts
load(['dat',sprintf('%4.4d',nt)])
phit(nt-st+1)=mean(phi(ix,2:ny0+1,iz));
end
phit=phit-mean(phit);  
cr(:)=ifft(phit);
phif(:)=abs(cr(:)); 

a(1)=phif(1);
a(2:Nt/2)=2*phif(2:Nt/2);

fid=['phif_LCO',sprintf('%4.4d',i)];
save(fid,'a')

plot(vf(1:Nt/2),a(1:Nt/2),'m-','Linewidth',1);





st=601;         %% H
nts=700;
Nt=nts-st+1;
at=ntp*tau*Nt;
nt02=Nt/2;
pi2=2.*pi;
pi2t=pi2/at; 
vf=zeros(Nt,1);
      for j=1:nt02+1
         vf(j)=pi2t*(j-1);
      end 
      for j=nt02+2:Nt
         jj=j-Nt;
         vf(j)=pi2t*(jj-1);
      end 
cr=zeros(Nt,1);
phif=zeros(Nt,1);
phit=phif;

for nt=st:nts
load(['dat',sprintf('%4.4d',nt)])
phit(nt-st+1)=mean(phi(ix,2:ny0+1,iz));
end
phit=phit-mean(phit);  
cr(:)=ifft(phit);
phif(:)=abs(cr(:)); 

a(1)=phif(1);
a(2:Nt/2)=2*phif(2:Nt/2);

fid=['phif_H',sprintf('%4.4d',i)];
save(fid,'a')

plot(vf(1:Nt/2),a(1:Nt/2),'r-','Linewidth',1);



st=741;         %% H-L
nts=840;
Nt=nts-st+1;
at=ntp*tau*Nt;
nt02=Nt/2;
pi2=2.*pi;
pi2t=pi2/at; 
vf=zeros(Nt,1);
      for j=1:nt02+1
         vf(j)=pi2t*(j-1);
      end 
      for j=nt02+2:Nt
         jj=j-Nt;
         vf(j)=pi2t*(jj-1);
      end 
cr=zeros(Nt,1);
phif=zeros(Nt,1);
phit=phif;

for nt=st:nts
load(['dat',sprintf('%4.4d',nt)])
phit(nt-st+1)=mean(phi(ix,2:ny0+1,iz));
end
phit=phit-mean(phit);  
cr(:)=ifft(phit);
phif(:)=abs(cr(:)); 

a(1)=phif(1);
a(2:Nt/2)=2*phif(2:Nt/2);

fid=['phif_HL',sprintf('%4.4d',i)];
save(fid,'a')

plot(vf(1:Nt/2),a(1:Nt/2),'g-','Linewidth',1);



axis([0 2.5 0 0.5]);
xlabel('f c_s/a'); legend('L','L-H','H','H-L');
print(gcf,'-depsc','tildephi_f')

