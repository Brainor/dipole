%% 图2
iz=2;%2d
y=-dy:dy:aly;
plot_variables={'phi','pei','deni';...%变量名称
    'wi','G','N';...
    'vex','fluxG','fluxN';...
    'vey','fluxp','fluxn'}';
plot_titles={'\Phi','\widetilde{G}','\widetilde{N}','w','G','N','v_x','\widetilde{G}v_x','\widetilde{N}v_x','v_y','\widetilde{p}_e v_x','\widetilde{n}v_x'};
plot_variables={'phi','G';...%variables name
    'N','deni'}';%遍历顺序为先列后行
plot_titles={'\Phi','G','N','\widetilde{N}'};
subplot_grid=num2cell(size(plot_variables'));
nt_s=[33,157];
close(gcf);
figure;
nCol=4;nRow=2;
rowH=0.7/nRow;colW=0.7/nCol;
colX=0.05+linspace(0,0.96,nCol+1);colX=colX(1:end-1);colX(3:end)=colX(3:end)+0.02;
rowY=0.1+linspace(0.9,0,nRow+1);rowY=rowY(2:end);
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
%% 图1,extended
%% with N evolution, scan eta, omegapd
% Gamma_N=real(i*m*Phi*N')
% Gamma_p=real(i*m*Phi*p*)
omegapds=linspace(1,2,500);omegapds=linspace(1.2,2,1000/2);
etas=linspace(0.4,2,200);etas=linspace(0,2,400/2);
[X,Y]=meshgrid(etas,omegapds);
rhom=1/25;
gamma=5/3;
m=3;
m_perp=m;

gr=nan(size(X));
vectors=nan(size(X,1),size(X,2),3);%N, p, Phi
for i=1:length(omegapds)
    for j=1:length(etas)
        omegapd=omegapds(i);eta=etas(j);
        A=[0,1,omegapd*1/(1+eta)-1;
            -gamma,2*gamma,omegapd-gamma;
            0,-1/0.6/rhom^2/m_perp^2,0];
        [V,w]=eigs(A,1,'largestimag');%分别为特征向量和特征值
        if imag(w)>0
            gr(i,j)=m*max(w);
            vectors(i,j,:)=V;
        end
    end
end
gammaN=real(1i*m*vectors(:,:,3).*conj(vectors(:,:,1)));
gammap=real(1i*m*vectors(:,:,3).*conj(vectors(:,:,2)));
angleN=angle(vectors(:,:,3)./vectors(:,:,1));
anglep=angle(vectors(:,:,3)./vectors(:,:,2));
%% eta, omegapd, {gamma,omega}
h=subplot(2,2,1);
contour(X,Y,imag(gr),'LevelStep',1);
h.XLabel.String='\eta';
h.YLabel.Interpreter='latex';h.YLabel.String='$$\frac{\omega^\star_p}{\omega_d}$$';
h.Title.String='Growth Rate';colorbar;
[col,row]=find(isnan(gr'),1,'last');
annotation('textarrow',[0.24 0.17],[0.81 0.79],'String',{sprintf('$$(%.3g,%.3g)$$',etas(col),omegapds(row))},'Interpreter','latex','HeadStyle','vback3','HorizontalAlignment','left');
h=subplot(2,2,2);
contour(X,Y,real(gr),'LevelStep',1);
h.XLabel.String='\eta';
h.YLabel.Interpreter='latex';h.YLabel.String='$$\frac{\omega^\star_p}{\omega_d}$$';
h.Title.String='Frequency';colorbar;
omegapd_index=find(omegapds>=5/3,1);
subplot(2,2,3);
baseline=2/3;
baseline_text='$\frac{2}{3}$';
h=draw_plot({etas,[imag(gr(omegapd_index,:));real(gr(omegapd_index,:))]},'$$\frac{\omega^\star_p}{\omega_d}=1.67$$','\eta','\omega/\omega_d','Box','on');
line(h,'XData',baseline*ones(size(ix)),'YData',h.YLim,'LineStyle','--');
[h.XTick,sortIndex]=sort([baseline,h.XTick]);h.XTickLabel{find(sortIndex==1,1)}=baseline_text;h.TickLabelInterpreter='latex';
legend('Growth rate','Frequency','Box','off','Location','best')

%% m, omegapd, {gamma,omega}
omegapds=[2,1.9,1.8,1.7,1.6,1.5];%取特定omegapd
iy=zeros(size(omegapds));
etas=[0.67,1.5];
ms=1:10;
rhom=1/25;
gamma=5/3;
m=1;

gr=inf(length(omegapds),length(ms),length(etas));
for i=1:length(omegapds)
    for j=1:length(ms)
        for k=1:length(etas)
            omegapd=omegapds(i);m=ms(j);eta=etas(k);
            A=[0,1,omegapd*1/(1+eta)-1;
                -gamma,2*gamma,omegapd-gamma;
                0,-1/0.6/rhom^2/m^2,0];
            w=eigs(A,1,'largestimag');
            if imag(w)>0
                gr(i,j,k)=m*max(w);
            end
        end
    end
end
for i=1:length(etas)
    subplot(2,2,4)
    draw_plot({ms,imag(gr(:,:,i)),'.-','Linewidth',1},['$$\eta=',num2str(etas(i)),'$$'],'m','\gamma/\omega_d','Box','on');
end
lgd=legend(cellfun(@num2str,num2cell(omegapds),'UniformOutput',false),'Location','best','NumColumns',2,'Box','off');
lgd.Title.String='$$\frac{\omega_p^\star}{\omega_d}$$';lgd.Title.Interpreter='latex';lgd.Title.NodeChildren.Position=[-.1,0.3,0];
print('paper\fig_eigen_extended','-dpng');

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
if(sp)
    data2.pe0=[data2.pe0';s_p/tau];
    data2.p=[data2.p';s_p/tau.*x.^(4*gamma)];
end
if(sn)
    data2.den0=[data2.den0';s_den/tau];
    data2.n=[data2.n';2*s_den/tau.*x.^4];
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
ix=[.7,0.9,1,1.1];
ix_index=arrayfun(@(ix)find(x>ix,1),ix);
n=zeros(length(ix),nts);
for nt=1:nts
    data=load(sprintf('data/dat%4.4d.mat',nt),'deni');
    n(:,nt)=mean(data.deni(ix_index,:,2),2);
end
close;
subplot(1,2,1);
draw_plot({1:nts,n},'','t','N','Box','on');
subplot(1,2,2);
h=draw_plot({1:nts,(n-1).*repmat(ix',[1,nts]).^4},'','t','n','Box','on');

legend(cellstr(num2str(ix','$$x=%g$$')),'Interpreter','latex','Box','off','Location','best');
print('paper/fig_N_t','-dpng')
%% 图5
%% 随时间演化的flux
ix=[find(x>0.7,1),find(x>.9,1)];%最靠近x=1的位置
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
draw_plot({(nt_start:nts),[flux_pt(i,:);flux_nt(i,:)]},'','t','','Box','on');
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
    line(h,'XData',h.XLim,'YData',baseline{i}*ones(size(ix)),'LineStyle','--')
    [h.YTick,sortIndex]=sort([baseline{i},h.YTick]);h.YTickLabel{find(sortIndex==1,1)}=baseline_text{i};h.TickLabelInterpreter='latex';
    legend(cellstr(num2str(ix','$$x=%g$$')),'Interpreter','latex','Box','off');
end
print('paper/fig_eta','-dpng');
%% 图7
%% 非线性时间平均
ix1=find(x>1,1);%最靠近x=1的位置
ix2=find(x>0.7,1);
iz=2;
nt=36;
data=load(sprintf('data/dat%4.4d.mat',nt),'phi','pei','deni');
phim1=amplitude(data.phi(ix1,2:end-1,iz));
phim2=amplitude(data.phi(ix2,2:end-1,iz));
denm1=amplitude(data.deni(ix1,2:end-1,iz));
pem1=amplitude(data.pei(ix1,2:end-1,iz));
denm2=amplitude(data.deni(ix2,2:end-1,iz));
pem2=amplitude(data.pei(ix2,2:end-1,iz));

denm1(1)=inf;
denm2(1)=inf;
subplot(1,2,1)
draw_plot({0:59,[phim1(1:60);denm1(1:60);phim2(1:60);denm2(1:60)]},'','m','','YScale','log','Box','on');

phim1=zeros(1,floor(ny/2));
pem1=phim1;denm1=phim1;
phim2=phim1;
pem2=phim1;denm2=phim1;
for nt=150:nts
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
denm1(1)=inf;
denm2(1)=inf;
subplot(1,2,2)
draw_plot({0:59,[phim1(1:60);denm1(1:60);phim2(1:60);denm2(1:60)]},'','m','','YScale','log','Box','on');
lgd=legend('$$\Phi(x=1)$$','$$N(x=1)$$','$$\Phi(x=0.7)$$','$$N(x=0.7)$$');lgd.Interpreter='latex';lgd.Box='off';

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