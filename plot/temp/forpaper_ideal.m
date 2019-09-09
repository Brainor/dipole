%% 图2
iz=2;%2d
y=-dy:dy:aly;
plot_variables={'phi','N';...%variables name
    'G','pei'}';%遍历顺序为先列后行
plot_titles={'\Phi','N','G','\widetilde{G}'};
nt_s=[30,157];
close(gcf);
for i=1:length(nt_s)
    data=load(sprintf('data/dat%4.4d.mat',nt_s(i)),'phi','vex','deni','pei','wi','vey');
    if (isdeltaf)%设定作图的变量
        data.N=data.deni+den0;
        data.G=data.pei+pe0;
    else
        data.N=data.deni;data.deni=delt(data.deni);
        data.G=data.pei;data.pei=delt(data.pei);
    end
    data.p=data.G.*repmat(x',[1,size(data.deni,2),size(data.deni,3)]).^(4*gamma);
    data.n=data.N.*repmat(x',[1,size(data.deni,2),size(data.deni,3)]).^4;
    data.fluxG=data.vex.*data.pei;
    data.fluxN=data.vex.*data.deni;
    data.fluxp=data.vex.*delt(data.p,2);
    data.fluxn=data.vex.*delt(data.n,2);
    for j=1:numel(plot_variables)
        %         h=axes('Position',[colX(-mod(j,2)+2*i),rowY(ceil(j/2)),colW,rowH]);
        %         draw_pcolor(h,{data.(plot_variables{j}),y,x},['$$',plot_titles{j},'$$'],'y','x');
        draw_pcolor({4,2,2*j-2+i},{data.(plot_variables{j}),y,x},['$$',plot_titles{j},'$$'],'y','x');
        %         h.XLabel.FontSize=14;h.YLabel.FontSize=16;h.Title.FontSize=20;
    end
    %     h=axes('Position',[i/2-0.5,0.95,0.5,0.05],'Color','none','XColor','none','YColor','none');
    %     text(h,0.5,0,num2str(nt_s(i)),'FontSize',24,'HorizontalAlignment','center');
end
print('paper/fig_ivp','-dpng');
%% 图1 global eigen
isdiffusive=1;%考虑耗散项
D=.06;
ms=1:10;growthrate=zeros(1,length(ms));phim_eigen=[];
h_profile=@(h0,c,x)h0+(1-h0)/(xmax-1)^(c*(xmax-1)/(1-xmin))/(1-xmin)^c*(xmax-x).^(c*(xmax-1)/(1-xmin)).*(x-xmin).^c;
for i=1:length(ms)
    m=ms(i);
    %默认参数
    gamma=5/3;
    L0=0.9;
    xmin=L0/2.5; xmax=L0/0.6;
    nx0=500*2*2*2;
    rhoL=2;%rho_s/L
    rhoL=.04;
    
    alx=xmax-xmin;
    dx_eigen=alx/(nx0+1);
    x_eigen=xmin:dx_eigen:xmax;
    omegad=4./x_eigen(2:end-1);
    
    %平衡剖面设置
    
    den0=h_profile(0.1,0.5,x_eigen);%den0=ones(size(x));
    pe0=h_profile(0,1,x_eigen);
    
    den1=den0(2:end)-den0(1:end-1);%hN(0.5):dx:hN(N+0.5)
    pe1=pe0(2:end)-pe0(1:end-1);%hG(0.5):dx:hG(N+0.5)
    x1_eigen=xmin+dx_eigen/2:dx_eigen:xmax-dx_eigen/2;%x(0.5):dx:x(N+0.5)
    den1=h_profile(0.1,0.5,x1_eigen);%den1=ones(size(x1));
    pe1=h_profile(0,1,x1_eigen);
    
    %先记录ideal
    M=zeros(2*nx0);
    N=zeros(2*nx0);
    
    M(1:nx0,1:nx0)=eye(nx0);
    
    ai=0.7/dx_eigen^2.*x1_eigen(2:end-1).^-2.*den1(2:end-1);
    bi=-0.7/dx_eigen^2.*(x1_eigen(2:end).^-2.*den1(2:end)+x1_eigen(1:end-1).^-2.*den1(1:end-1))-0.6*m^2*den0(2:end-1).*x_eigen(2:end-1).^-4;
    ci=0.7/dx_eigen^2.*x1_eigen(2:end-1).^-2.*den1(2:end-1);
    M(nx0+1:2*nx0,nx0+1:2*nx0)=diag(bi)+diag(ai,-1)+diag(ci,1);
    
    N(1:nx0,nx0+1:2*nx0)=m/dx_eigen*diag((pe1(2:end)-pe1(1:end-1)).*x_eigen(2:end-1).^(4*gamma));
    % N(nx0+1:2*nx0,1:nx0)=m/rhoL^2*diag(omegad.*x(2:end-1).^(-4*gamma));
    N(nx0+1:2*nx0,1:nx0)=m/rhoL^2*diag(omegad.*x_eigen(2:end-1).^(-4));
    
    M=sparse(M);
    N=sparse(N);
    
    [V,w]=eigs(N,M,3,'largestimag');
    if m==6
        subplot(2,2,1)
        draw_plot({1./x_eigen(2:end-1),real(V(nx0+1:2*nx0,:))},['$$m=',num2str(m),'$$, $$\rho_\star=',num2str(rhoL),'$$'],'x','\Re\{\Phi\}','XAxisLocation','origin');
        lgd=legend(cellfun(@num2str,num2cell(imag(diag(w))),'UniformOutput',false),'Location','best','Box','off');
        lgd.Title.String='$$\gamma$$';lgd.Title.Interpreter='latex';
        subplot(2,2,2)
        draw_plot({x_eigen(2:end-1),[den0(2:end-1);pe0(2:end-1)]},'Equilibruim Profile','x');legend('$$h_N$$','$$h_G$$','Interpreter','latex','Box','off','location','best');
        h=subplot(2,2,3);
        X=x_eigen(2:end-1)'*cos(linspace(0,2*pi,100));
        Y=x_eigen(2:end-1)'*sin(linspace(0,2*pi,100));
        phi=real(V(nx0+1:2*nx0,1));
        C=phi*cos(m*linspace(0,2*pi,100));
        draw_pcolor(h,{C,X,Y});
    end
    if isdiffusive
        growthrate(i)=imag(w(1))-D*m^2;
    else
        growthrate(i)=imag(w(1));
    end
    phim_eigen(i,:)=real(V(nx0+1:end,1));
end
subplot(2,2,4)
draw_plot({ms,growthrate},{'growth rate',sprintf('$$D=%g$$',D)},'m','\gamma');
print('paper/fig_eigen_ideal','-dpng');
%% 图3
%% 环向和时间平均profile
nt_s=150:200;
plot_variables={'pe0','den0';...%variables name
    'fluxp','fluxn';...
    'p','n'}';%遍历顺序为先列后行
plot_titles={'G','N','\Gamma_G','\Gamma_N','p_e','n'};
subplot_grid=num2cell(size(plot_variables'));
iz=2;
data2=struct;
for nt=nt_s
    data=load(sprintf('data/dat%4.4d.mat',nt),'phi','pei','vex','deni','wi');
    if (isdeltaf)
        data.pe0=mean(pe0(:,:,iz)-bc_p,2);
        data.den0=mean(den0(:,:,iz)-bc_n,2);
    else
        data.pe0=mean(data.pei(:,:,iz)-bc_p,2);
        data.den0=mean(data.deni(:,:,iz)-bc_n,2);
    end
    data.phi=mean(data.phi(:,:,iz),2);data.wi=mean(data.wi(:,:,iz),2);
    data.fluxp=mean(data.vex(:,:,iz).*data.pei(:,:,iz),2);
    data.fluxn=mean(data.vex(:,:,iz).*data.deni(:,:,iz),2);
    data.n=mean(data.den0.*repmat(x',[1,size(data.deni,2)]).^4,2);
    data.p=mean(data.pe0.*repmat(x',[1,size(data.deni,2)]).^(4*gamma),2);
    for i=1:numel(plot_variables)
        data2.(plot_variables{i})(:,nt-nt_s(1)+1)=data.(plot_variables{i});
    end
end
for i=1:numel(plot_variables)%做时间平均
    data2.(plot_variables{i})=mean(data2.(plot_variables{i}),2);
end
sp=1;sn=1;
if(sp)
    data2.pe0=[data2.pe0';s_p/tau/5];
    data2.p=[data2.p';s_p/tau.*x.^(4*gamma)/5];
end
if(sn)
    data2.den0=[data2.den0';s_den/tau/6];
    data2.n=[data2.n';s_den/tau.*x.^4];
end%加上源
% 横坐标为x
for i=1:numel(plot_variables)
    subplot(subplot_grid{:},i);
    h=draw_plot({x,data2.(plot_variables{i})},['$$',plot_titles{i},'$$'],'x','','Box','on');
    if(length(h.Children)==2),legend(h.Children(1),'source','Box','off','Location','best');end
end
print('paper/fig_ivp_flux','-dpng');
%% 图4
%% 环向平均N随时间演化, 给定x
ix=[0.7,.8,.9,1];
ix_index=arrayfun(@(ix)find(x>ix,1),ix);
n=zeros(length(ix),nts);G=n;
for nt=1:nts
    data=load(sprintf('data/dat%4.4d.mat',nt),'deni','pei');
    n(:,nt)=mean(data.deni(ix_index,:,2),2);
    G(:,nt)=mean(data.pei(ix_index,:,2),2);
end
subplot(2,2,1)
draw_plot({1:nts,n},'','t','N','Box','on');
subplot(2,2,2)
draw_plot({1:nts,(n-1).*repmat(ix',[1,nts]).^4},'','t','n','Box','on');
legend(cellstr(num2str(ix','$$x=%g$$')),'Interpreter','latex','Box','off','Location','best');
subplot(2,2,3)
draw_plot({1:nts,G},'','t','G','Box','on');
subplot(2,2,4)
draw_plot({1:nts,(G-1).*repmat(ix',[1,nts]).^(4*gamma)},'','t','p_e','Box','on');
print('paper/fig_N_t','-dpng')
%% 图5
%% 随时间演化的flux
ix=[find(x>0.9,1),find(x>0.7,1)];%最靠近x=1的位置
iy=ny/2;
iz=2;
nt_start=1;
flux_pt=zeros(length(ix),nts-nt_start+1);
flux_nt=flux_pt;
for nt=nt_start:nts
    data=load(sprintf('data/dat%4.4d.mat',nt),'pei','vex','deni');
    %     flux_pt(:,nt+1-nt_start)=data.pei(ix,iy,iz).*data.vex(ix,iy,iz);
    %     flux_nt(:,nt+1-nt_start)=data.deni(ix,iy,iz).*data.vex(ix,iy,iz);
    flux_pt(:,nt+1-nt_start)=mean(data.pei(ix,:,iz).*data.vex(ix,:,iz),2);
    flux_nt(:,nt+1-nt_start)=mean(data.deni(ix,:,iz).*data.vex(ix,:,iz),2);
end
%%
close
for i=1:length(ix)
    subplot(1,length(ix),i);
    draw_plot({(nt_start:nts),[flux_pt(i,:);flux_nt(i,:)]},'','t','\Gamma','Box','on');
end
legend('$$\Gamma_G$$','$$\Gamma_N$$','Interpreter','latex','Box','off');
print('paper/fig_flux_t','-dpng');
%% 图6
%% eta和omegapd特定x位置随时间演化
ix=[0.7,0.9];
ix_index=arrayfun(@(ix)find(x>ix,1),ix);
baseline={5/3,2/3};
baseline_text={'$\frac{5}{3}$','$\frac{2}{3}$'};
plot_variables={'omegap','eta'}';%variables name, 遍历顺序为先列后行
plot_titles={'\frac{\omega_p^*}{\omega_d}','\eta'};
subplot_grid=num2cell(size(plot_variables'));
data2=struct;
for nt=1:nts
    data=load(sprintf('data/dat%4.4d.mat',nt));
    if isdeltaf
        pe0=pe0(:,:,2)';%%pe0为平均量
        den0=den0(:,:,2)';
    else
        pe0=mean(data.pei(:,:,2),2)';
        den0=mean(data.deni(:,:,2),2)';
    end
    dpe=gradient(pe0,dx);%dG/dx
    dden=gradient(den0,dx);%dN/dx
    eta=(dpe./pe0+4*gamma./x)./(dden./den0+4./x)-1;
    data2.eta(:,nt)=eta(ix_index);
    omegap=dpe./pe0/4.*x+gamma;
    data2.omegap(:,nt)=omegap(ix_index);
end
%%
close
for i=1:numel(plot_variables)
    subplot(subplot_grid{:},i);
    h=draw_plot({1:nts,data2.(plot_variables{i})},'','t',plot_titles{i},'Box','on');
    line(h,'XData',h.XLim,'YData',baseline{i}*ones(1,2),'LineStyle','--')
    [h.YTick,sortIndex]=sort([baseline{i},h.YTick]);h.YTickLabel{find(sortIndex==1,1)}=baseline_text{i};h.TickLabelInterpreter='latex';
    legend(cellstr(num2str(ix','$$x=%g$$')),'Interpreter','latex','Box','off');
end
print('paper/fig_eta','-dpng');
%% 图7
%% 非线性时间平均
ix1=find(x>1,1);%最靠近x=1的位置
ix2=find(x>0.7,1);
iz=2;
nt=30;
data=load(sprintf('data/dat%4.4d.mat',nt),'phi','pei','deni');
phim1=amplitude(data.phi(ix1,2:end-1,iz));
phim2=amplitude(data.phi(ix2,2:end-1,iz));
denm1=amplitude(data.deni(ix1,2:end-1,iz));
pem1=amplitude(data.pei(ix1,2:end-1,iz));
denm2=amplitude(data.deni(ix2,2:end-1,iz));
pem2=amplitude(data.pei(ix2,2:end-1,iz));

denm1(1)=inf;pem1(1)=inf;
denm2(1)=inf;pem2(1)=inf;
subplot(1,2,1)
draw_plot({0:59,[phim1(1:60);pem1(1:60);phim2(1:60);pem2(1:60)]},'','m','','YScale','log','Box','on');

phim1=zeros(1,floor(ny/2));
pem1=phim1;denm1=phim1;
phim2=phim1;
pem2=phim1;denm2=phim1;
for nt=50:nts
    data=load(sprintf('data/dat%4.4d.mat',nt));
    
    phim1=amplitude(data.phi(ix1,2:end-1,iz))+phim1;
    phim2=amplitude(data.phi(ix2,2:end-1,iz))+phim2;
    if (isdeltaf)
        pem=abs(fft(data.pei(ix,:,iz)+pe0(ix,:,iz)))/ny+pem;
        denm1=abs(fft(data.deni(ix,:,iz)+den0(ix,:,iz)))/ny+denm1;
    else
        pem1=amplitude(data.pei(ix1,2:end-1,iz))+pem1;
        denm1=amplitude(data.deni(ix1,2:end-1,iz))+denm1;
        pem2=amplitude(data.pei(ix2,2:end-1,iz))+pem2;
        denm2=amplitude(data.deni(ix2,2:end-1,iz))+denm2;
    end
end
phim1=phim1/(nts-150+1);
pem1=pem1/(nts-150+1);
denm1=denm1/(nts-150+1);
phim2=phim2/(nts-150+1);
pem2=pem2/(nts-150+1);
denm2=denm2/(nts-150+1);
denm1(1)=inf;pem1(1)=inf;
denm2(1)=inf;pem2(1)=inf;
subplot(1,2,2)
draw_plot({0:59,[phim1(1:60);pem1(1:60);phim2(1:60);pem2(1:60)]},'','m','','YScale','log','Box','on');
legend('$$\Phi(x=1)$$','$$G(x=1)$$','$$\Phi(x=0.7)$$','$$G(x=0.7)$$','Interpreter','latex','Box','off','Location','best');

print('paper/fig_spectrum','-dpng');

function vm=amplitude(v)
% AMPLITUDE calculate the amplitude of phase
% v is the real vector
% vm has length of floor(n/2)+1, with phase number 0:floor(n/2)+1
n=length(v);
vm=abs(fft(v))/n;
vm(2:floor((n+1)/2))=2*vm(2:floor((n+1)/2));
vm(floor(n/2)+2:end)=[];
end
