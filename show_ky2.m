global gamma
close all
m=[0:ny0/2,floor(-ny0/2)+1:-1];

     cr=zeros(ny0,1);
     cr2=zeros(ny0,1);
phik=zeros(ny0,nts);
denk=zeros(ny0,nts);

n2=4; %choose ky range.
i=30; %choose x location.

for nt=1:nts
    nt
load(['dat',sprintf('%4.4d',nt)])
%densave(:,:,nt)=deni(:,:,2);



cr(:)=ifft(phi(i,2:ny0+1,2)); 
%cr(:)=ifft(tei(i,2:ny0+1,2)); 
cr2(:)=ifft(pei(i,2:ny0+1,2)/x(i).^(-4*gamma)); 

phik(:,nt)=abs(cr(:));
denk(:,nt)=abs(cr2(:));

a(1)=phik(1,nt);
a(2:ny0/n2)=2*phik(2:ny0/n2,nt);
semilogy(m(1:ny0/n2),a,'-','Linewidth',1)
hold

a(1)=denk(1,nt);
a(2:ny0/n2)=2*denk(2:ny0/n2,nt);
semilogy(m(1:ny0/n2),a,'-r','Linewidth',1)
%semilogy(phik(1:ny0,nt),'-b','Linewidth',1)
xlabel('k_y \rho_s'); legend('\phi','p','Location','best');
print(gcf,'-dpng',sprintf('ky%4.4d',nt)) 
hold off

end
close all

%%

