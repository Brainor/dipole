% use this in linear part
% given x & z, mean in y, check the time varing
% std(\tilde{p(ix,:,iz)}_e) time varing 
cntnu=1;
petyt=[];
petyt2=[];
petyt3=[];
petyt4=[];
ix=13;iz=6;

lsdat=dir('dat*.mat');
ndat=length(lsdat);

for nt=1:ndat
    load(['dat',sprintf('%4.4d',nt)]);
    pey =pei(ix ,:,iz);    petyt (end+1)=std(pey, 1);
end


lt=length(petyt);
t=(1:lt)/10;

figure;
set(gca,'fontsize',14);
plot(t,petyt);
title('$\tilde{p}_e -t(13x)$','interpreter','latex');
xlabel('t');ylabel('$\tilde{p}_e$','interpreter','latex');
print(gcf,'-dpng',sprintf('pet2t%2.2dx',ix));
close

% log p -t
figure;
set(gca,'fontsize',14);
semilogy(t,petyt);
title('$\log\ \tilde{p}_e -t(13x)$','interpreter','latex');
xlabel('t');ylabel('$\tilde{p}_e$','interpreter','latex');
print(gcf,'-dpng',sprintf('pet2t%2.2dxlogy',ix));
close

% get the growth rate of pe tilde
Dt=tau*ntp;
gpe=firstlog(petyt,Dt);
figure;
plot(t,gpe)
set(gca,'fontsize',14)
xlabel('t')
ylabel('\gamma')
title('$\gamma of \tilde{p}_e -t$','interpreter','latex');
print(gcf,'-dpng',sprintf('pegamma%2.2dx',ix));
close
