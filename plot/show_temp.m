%% 约束时间
nt_s=150:nts;
N=zeros(2,length(nt_s));G=N;%第一行为数据, 第二行为约束时间
for i=1:length(nt_s)
    data=load(sprintf('data/dat%4.4d.mat',nt_s(i)),'deni','pei','deni');
    N(1,i)=sum(data.deni(2:end-1,2:end-1,2),'all')*dx*dy;G(1,i)=sum(data.pei(2:end-1,2:end-1,2),'all')*dx*dy;
end
N(2,:)=N(1,:)./(sum(source_den(2:end-1,2:end-1,2),'all')*ntp*dx*dy-gradient(N(1,:),tau*ntp));
G(2,:)=G(1,:)./(sum(source_p(2:end-1,2:end-1,2),'all')*ntp*dx*dy-gradient(G(1,:),tau*ntp));
h=subplot(1,2,1);
draw_plot({h,nt_s,[N(2,:);G(2,:)]},'$$\tau=\frac{E}{P_{in}-dE/dt}$$','t','\tau/\frac{L_0}{\rho_*c_s}');
legend({'$\tau_N$','$\tau_G$'},'Interpreter','latex','Box','off','Location','best');
%% diffusion loss rate
nt=150;
data=load(sprintf('data/dat%4.4d.mat',nt),'deni','pei','deni');
N=-dif*gradient(gradient(sum(data.deni(:,2:end-1,2),2),dx),dx);
G=-dif*gradient(gradient(sum(data.pei(:,2:end-1,2),2),dx),dx);
h=subplot(1,2,2);
draw_plot({h,x,[N,G]},'$$-D\frac{d^2N}{dx^2}$$','x');
legend({'$N$','$G$'},'Interpreter','latex','Box','off','Location','best');
%% 各种分量时间平均图
x_s=1.3;
nt_start=170;nts=170;
ix=arrayfun(@(x_s)find(x>x_s,1),x_s);
phim=zeros(length(x_s),ny0);pem=phim;denm=phim;vxm=phim;fluxnm=phim;
for nt=nt_start:nts
    data=load(sprintf('data/dat%4.4d.mat',nt),'phi','deni','pei','vex');
    for i=1:length(x_s)
        phim(i,:)=data.phi(ix(i),2:end-1,iz)+phim(i,:);
        denm(i,:)=data.deni(ix(i),2:end-1,iz)+denm(i,:);
        pem(i,:)=data.pei(ix(i),2:end-1,i)+pem(i,:);
        vxm(i,:)=data.vex(ix(i),2:end-1,iz)+vxm(i,:);
        fluxnm(i,:)=data.vex(ix(i),2:end-1,iz).*data.deni(ix(i),2:end-1,iz)+fluxnm(i,:);
    end
end
% phim=phim/(nts-nt_start+1);pem=pem/(nts-nt_start+1);denm=denm/(nts-nt_start+1);vxm=vxm/(nts-nt_start+1);fluxnm=fluxnm/(nts-nt_start+1);
phim=phim/max(abs(phim));pem=pem/max(abs(pem));denm=denm/max(abs(denm));vxm=vxm/max(abs(vxm));fluxnm=fluxnm/max(abs(fluxnm));
figure
h=draw_plot({0:dy:2*pi-dy,[phim;denm;vxm;fluxnm]},'','y','','Box','on','XAxisLocation','origin','FontSize',12,'OuterPosition',[0,0,1,1]);h.LooseInset=h.TightInset;
legend({'$\Phi$','$N$','$v_x$','$\Gamma_N$'},'Interpreter','latex','Box','off','Location','best');
%% 各分量时间演化
x_s=1;
y_s=pi;
ix=arrayfun(@(x_s)find(x>x_s,1),x_s);
iy=ny0/2;
nt_start=1;
phim=zeros(length(ix),length(nt_start:nts));pem=phim;denm=phim;vxm=phim;fluxnm=phim;
for nt=nt_start:nts
    data=load(sprintf('data/dat%4.4d.mat',nt),'phi','deni','pei','vex');
    for i=1:length(x_s)
        phim(i,nt-nt_start+1)=data.phi(ix(i),iy,iz);
        denm(i,nt-nt_start+1)=data.deni(ix(i),iy,iz);
        pem(i,nt-nt_start+1)=data.pei(ix(i),iy,iz);
        vxm(i,nt-nt_start+1)=data.vex(ix(i),iy,iz);
        fluxnm(i,nt-nt_start+1)=data.vex(ix(i),iy,iz).*data.deni(ix(i),iy,iz);
    end
end
figure;
draw_plot({nt_start:nts,[denm;pem]},'','t','','FontSize',12,'OuterPosition',[0,0,1,1]);h.LooseInset=h.TightInset;
legend({'$N$','$G$'},'Interpreter','latex','Box','off','Location','best');
%% 某个点的通量, 沿着y方向
x_s=1;
ix=arrayfun(@(x_s)find(x>x_s,1),x_s);
nt=170;
data=load(sprintf('data/dat%4.4d.mat',nt),'phi','deni','pei','vex');
data.fluxN=data.deni.*data.vex;data.fluxdN=delt(data.deni,2).*data.vex;
h=subplot(1,2,1);
draw_plot({h,0:dy:2*pi-dy,[data.deni(ix,2:end-1,2);data.pei(ix,2:end-1,2)]},'','y','','FontSize',12,'OuterPosition',[0,0,1/2,1]);h.LooseInset=h.TightInset;
legend({'$N$','$G$'},'Interpreter','latex','Box','off','Location','best');
h=subplot(1,2,2);
draw_plot({h,0:dy:2*pi-dy,[data.fluxN(ix,2:end-1,2);data.fluxdN(ix,2:end-1,2)]},'','y','','FontSize',12,'OuterPosition',[1/2,0,1/2,1]);h.LooseInset=h.TightInset;
legend({'$\Gamma_N$','$\Gamma_{\delta N}$'},'Interpreter','latex','Box','off','Location','best');
annotation('textbox',[h.Position(1),h.Position(4),0.1,0.1],'String',num2str(x_s(1),'$$x=%g$$'),'FitBoxToText','on','Interpreter','latex','EdgeColor','none','FontSize',12);