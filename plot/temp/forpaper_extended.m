set(groot,'defaultAxesFontSize',14)
set(groot,'defaultAxesFontSize','remove')
%% 图1,extended
%% with N evolution, scan eta, omegapd
% Gamma_N=real(i*m*Phi*N')
% Gamma_p=real(i*m*Phi*p*)
omegapds=linspace(1,2,500);omegapds=linspace(1.3,2,1000/2);
etas=linspace(0.4,2,200);etas=linspace(0.2,1.6,400/2);
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
close
h_s = tight_subplot(1,2,[.03 .01],[.1 .03],[.09 .01]);
% h=subplot(1,2,1);
contour(h_s(1),X,Y,imag(gr),'LevelStep',1);h_s(1).FontSize=12;%h.OuterPosition=[(1-1)/2,0,1/2,1];
h_s(1).XLabel.String='\eta';h_s(1).XLabel.Position(2)=h_s(1).XLabel.Position(2)+0.02;
h_s(1).YLabel.Interpreter='latex';h_s(1).YLabel.String='$$\omega^\star_p/\omega_d$$';h_s(1).YLabel.Position(1)=h_s(1).YLabel.Position(1)+0.01;
colorbar(h_s(1));%h.LooseInset=h.TightInset;
% h=subplot(1,2,2);
% contour(h_s(2),X,Y,real(gr),'LevelStep',1);h_s(2).FontSize=12;%h.OuterPosition=[(2-1)/2,0,1/2,1];
contour(h_s(2),X,Y,real(gr),'LevelStep',1,'LevelList',[-6:-1,1:12]);h_s(2).FontSize=12;%h.OuterPosition=[(2-1)/2,0,1/2,1];
hold on
contour(h_s(2),X,Y,real(gr),'LevelStep',1,'LevelList',0,'LineStyle','--');h_s(2).FontSize=12;%h.OuterPosition=[(2-1)/2,0,1/2,1];
hold off
h_s(2).XLabel.String='\eta';h_s(2).YTickLabel='';h_s(2).XLabel.Position(2)=h_s(2).XLabel.Position(2)+0.02;
colorbar(h_s(2));%h.LooseInset=h.TightInset;
print('paper\fig_eigen_extended','-dpng');
print('paper\fig_eigen_extendeds','-depsc','-tiff');
%% 图2
iz=2;%2d
plot_variables={'phi','pei','deni';...%变量名称
    'wi','G','N';...
    'vex','fluxG','fluxN';...
    'vey','fluxp','fluxn'}';
plot_titles={'\Phi','\widetilde{G}','\widetilde{N}','w','G','N','v_x','\widetilde{G}v_x','\widetilde{N}v_x','v_y','\widetilde{p}_e v_x','\widetilde{n}v_x'};
plot_variables={'phi','G','N','vex'}';%遍历顺序为先列后行
plot_titles={'\Phi','G','N','v_x'};
subplot_grid=num2cell(size(plot_variables'));
nt_s=[66,155];
nt_s=[84,142];nt_s=[102,198];
figure;
nCol=4;nRow=2;
rowH=0.7/nRow;colW=0.7/nCol;
colX=0.05+linspace(0,0.96,nCol+1);colX=colX(1:end-1);colX(3:end)=colX(3:end)+0.02;
rowY=0.1+linspace(0.9,0,nRow+1);rowY=rowY(2:end);
h_s = tight_subplot(4,2,[.05 .015],[.09 .03],[.07 0]);
%%
for i=2:length(nt_s)
    data=load(sprintf('data/dat%4.4d.mat',nt_s(i)),'phi','vex','deni','pei','wi','vey');
    if (isdeltaf)%设定作图的变量
        data.N=data.deni+den0;
        data.G=data.pei+pe0;
    else
        data.N=data.deni-bc_n;data.deni=delt(data.deni);
        data.G=data.pei-bc_p;data.pei=delt(data.pei);
    end
    data.p=data.G.*repmat(x',[1,size(data.deni,2),size(data.deni,3)]).^(4*gamma);
    data.n=data.N.*repmat(x',[1,size(data.deni,2),size(data.deni,3)]).^4;
    data.fluxG=data.vex.*data.pei;
    data.fluxN=data.vex.*data.deni;
    data.fluxp=data.vex.*delt(data.p,2);
    data.fluxn=data.vex.*delt(data.n,2);
    
    data.dG=delt(data.pei);
    for j=1:numel(plot_variables)
        %         h=axes('Position',[colX(-mod(j,2)+2*i),rowY(ceil(j/2)),colW,rowH]);
        %         draw_pcolor(h,{data.(plot_variables{j}),y,x},['$$',plot_titles{j},'$$'],'y','x');
        
        %         h=draw_pcolor({4,2,2*j-2+i},{data.(plot_variables{j}),y,x},['$$',plot_titles{j},'$$'],'','x','FontSize',12);
        %         if j==4,h.XLabel.String='$$y$$';h.XLabel.Interpreter='latex';end
        %         else,h.XTick=[];end
        %         h.XLabel.FontSize=14;h.YLabel.FontSize=16;h.Title.FontSize=20;
        h=draw_pcolor(h_s(2*j-2+i),{data.(plot_variables{j}),-dy:dy:aly,x},['$$',plot_titles{j},'$$'],'','','FontSize',12);
        h.XLabel.Position(2)=h.XLabel.Position(2)+0.15;h.Title.VerticalAlignment='bottom';
        h.Title.Position(2)=h.Title.Position(2)-0.13;h.YLabel.Position(1)=h.YLabel.Position(1)+0.4;
        if(j==4),h.Title.Position(2)=h.Title.Position(2)+0.04;end
%         h.Title.Position
    end
    %     h=axes('Position',[i/2-0.5,0.95,0.5,0.05],'Color','none','XColor','none','YColor','none');
    %     text(h,0.5,0,num2str(nt_s(i)),'FontSize',24,'HorizontalAlignment','center');
end
set(h_s(1:6),'XTickLabel','');set(h_s(2:2:8),'YTickLabel','');
for i=1:2:8,h_s(i).YLabel.String='$$x$$';h_s(i).YLabel.Interpreter='latex';end
for i=[7,8],h_s(i).XLabel.String='$y$';h_s(i).XLabel.Interpreter='latex';end
%%
print('paper/fig_ivp_2','-dpng');
print('paper/fig_ivp','-depsc','-tiff');

%% 图3
%% 环向和时间平均profile
nt_s=198:200;nt_s=180:200;
plot_variables={'pe0','den0';...%variables name
    'fluxp','fluxn';...
    'p','n'}';%遍历顺序为先列后行
plot_titles={'G','N','\Gamma_G','\Gamma_N','p_e','n'};
plot_variables={'pe0','den0';...%variables name
    'p','n'}';%遍历顺序为先列后行
plot_titles={'G','N','p','n'};
subplot_grid=num2cell(size(plot_variables'));
iz=2;
data2=struct;
for nt=nt_s
    data=load(sprintf('data/dat%4.4d.mat',nt),'phi','pei','vex','deni','wi');
    if (isdeltaf)
        data.pe0=mean(pe0(:,2:end-1,iz),2);
        data.den0=mean(den0(:,2:end-1,iz),2);
    else
        data.pe0=mean(data.pei(:,2:end-1,iz),2);
        data.den0=mean(data.deni(:,2:end-1,iz),2);
    end
    data.phi=mean(data.phi(:,2:end-1,iz),2);data.wi=mean(data.wi(:,2:end-1,iz),2);
    data.fluxp=mean(data.vex(:,2:end-1,iz).*data.pei(:,2:end-1,iz),2);
    data.fluxn=mean(data.vex(:,2:end-1,iz).*data.deni(:,2:end-1,iz),2);
    data.n=mean((data.den0-bc_p(1)).*repmat(x',[1,size(data.deni,2)]).^4,2);
    data.p=mean((data.pe0-bc_n(1)).*repmat(x',[1,size(data.deni,2)]).^(4*gamma),2);
    for i=1:numel(plot_variables)
        data2.(plot_variables{i})(:,nt-nt_s(1)+1)=data.(plot_variables{i});
    end
end
for i=1:numel(plot_variables)%做时间平均
    data2.(plot_variables{i})=mean(data2.(plot_variables{i}),2);
end
if(sp)
    data2.pe0=[data2.pe0';s_p/tau/2+bc_p(1)];
    %     data2.p=[data2.p';s_p/tau.*x.^(4*gamma)];
end
if(sn)
    data2.den0=[data2.den0';s_den/tau/2+bc_n(1)];
    %     data2.n=[data2.n';2*s_den/tau.*x.^4/2];
end%加上源
%%
% 横坐标为x
close
% h_s = tight_subplot(2,2,[.05 .08],[.07 .05],[.06 .02]);
for i=1:numel(plot_variables)
    subplot(subplot_grid{:},i);
    h=draw_plot({x,data2.(plot_variables{i})},'','',plot_titles{i},'Box','on','FontSize',12);
    if(length(h.Children)==2),legend(h.Children(1),'Source','Box','off','Location','best');end
    if i>2,h.XLabel.String='$$x$$';h.XLabel.Interpreter='latex';end
    annotation('textbox',[h.Position(1),min(h.Position(2)+h.Position(4)+0.02,1),0,0],'String',['(',char(96+i),')'],'FitBoxToText','on','EdgeColor','none','FontSize',12);
    %     h=draw_plot({h_s(i),x,data2.(plot_variables{i})},['$$',plot_titles{i},'$$'],'','','Box','on','FontSize',12);
    %     if(length(h.Children)==2),legend(h.Children(1),'Source','Box','off','Location','best');end
    %     h.Title.Units='normalized';h.Title.Position(2)=h.Title.Position(2)-0.05;
end
% set(h_s(1:2),'XTickLabel','');
% for i=3:4
%     h_s(i).XLabel.String='$$x$$';h_s(i).XLabel.Interpreter='latex';
%     h_s(i).XLabel.Units='normalized';h_s(i).XLabel.Position(2)=h_s(i).XLabel.Position(2)+0.1;
% end
print('paper/fig_ivp_flux','-dpng');
print('paper/fig_ivp_flux','-depsc');
%% 图4图5合并
%% 图4
%% 随时间演化的flux
x_s=[0.8,0.9];
ix=arrayfun(@(x_s)find(x>x_s,1),x_s);
iy=ny/2;
iz=2;
nt_start=1;
flux_pt=zeros(length(ix),nts-nt_start+1);
flux_nt=flux_pt;
for nt=nt_start:nts
    data=load(sprintf('data/dat%4.4d.mat',nt),'pei','vex','deni');
    %     flux_pt(:,nt+1-nt_start)=data.pei(ix,iy,iz).*data.vex(ix,iy,iz);
    %     flux_nt(:,nt+1-nt_start)=data.deni(ix,iy,iz).*data.vex(ix,iy,iz);
    flux_pt(:,nt+1-nt_start)=mean(data.pei(ix,2:end-1,iz).*data.vex(ix,2:end-1,iz),2);
    flux_nt(:,nt+1-nt_start)=mean(data.deni(ix,2:end-1,iz).*data.vex(ix,2:end-1,iz),2);
end
%% 图5
%% 环向平均N随时间演化, 给定x
ix=[0.8,0.9,1];
ix_index=arrayfun(@(ix)find(x>ix,1),ix);
n=zeros(length(ix),nts);g=n;
for nt=1:nts
    data=load(sprintf('data/dat%4.4d.mat',nt),'deni','pei');
    n(:,nt)=mean(data.deni(ix_index,2:end-1,2),2);
    g(:,nt)=mean(data.pei(ix_index,2:end-1,2),2);
end
%%
close
for i=1:2
    h=subplot(2,2,i);
    draw_plot({(nt_start:nts)*tau*ntp,[flux_nt(i,:);flux_pt(i,:)]},'','','','Box','on','FontSize',12);%,'OuterPosition',[(i-1)/length(ix),0.5,1/length(ix),.5]);h.LooseInset=h.TightInset;
    h.Children(1).LineStyle='--';h.YLim(2)=0.1;h.XLim(1)=0;h.XTick=0:.2:1;
    legend('$$\Gamma_N$$','$$\Gamma_G$$','Interpreter','latex','Box','off','Location','best');
    annotation('textbox',[h.Position(1),h.Position(2),0.1,0.1],'String',['$$x=',num2str(x_s(i)),'$$'],'FitBoxToText','on','Interpreter','latex','EdgeColor','none','FontSize',12);
    annotation('textbox',[h.Position(1),h.Position(2)+h.Position(4)+0.02,0,0],'String',['(',char(96+i),')'],'FitBoxToText','on','EdgeColor','none','FontSize',12);
end
h=subplot(2,2,3);
draw_plot({(1:nts)*tau*ntp,n},'','t','N','Box','on','FontSize',12,'XTick',0:.2:1);%,'OuterPosition',[0,0,.5,.5]);h.LooseInset=h.TightInset;
% h.Title.Units='normalized';h.Title.Position(2)=h.Title.Position(2)-0.02;
h.Children(1).Color=h.ColorOrder(4,:);h.XLim(1)=0;
annotation('textbox',[h.Position(1),h.Position(2)+h.Position(4)+0.02,0,0],'String',['(','c',')'],'FitBoxToText','on','EdgeColor','none','FontSize',12);
% legend(cellstr(num2str(ix','$$x=%g$$')),'Interpreter','latex','Box','off','Location','best');h.XLabel.String='t';h.XLabel.Interpreter='latex';
h=subplot(2,2,4);
% draw_plot({(1:nts)*tau*ntp,(n).*repmat(ix',[1,nts]).^4},'$$n$$','t','','Box','on','FontSize',12,'OuterPosition',[(2-1)/2,0,1/2,1]);h.LooseInset=h.TightInset;
draw_plot({(1:nts)*tau*ntp,g},'','t','G','Box','on','FontSize',12,'XTick',0:.2:1);%,'OuterPosition',[0.5,0,0.5,.5]);h.LooseInset=h.TightInset;
% h.Title.Units='normalized';h.Title.Position(2)=h.Title.Position(2)-0.02;
h.Children(1).Color=h.ColorOrder(4,:);h.XLim(1)=0;
legend(cellstr(num2str(ix','$$x=%g$$')),'Interpreter','latex','Box','off','Location','best');%h.XLabel.String='t';h.XLabel.Interpreter='latex';
annotation('textbox',[h.Position(1),h.Position(2)+h.Position(4)+0.02,0,0],'String',['(','d',')'],'FitBoxToText','on','EdgeColor','none','FontSize',12);
%% tight plot
close
h_s = tight_subplot(2,2,[.035 .08],[.09 0.01],[.1 0.01]);
for i=1:2
    h=draw_plot({h_s(i),(nt_start:nts)*tau*ntp,[flux_nt(i,:);flux_pt(i,:)]},'','','','Box','on','FontSize',12,'XTickLabel','');%,'OuterPosition',[(i-1)/length(ix),0.5,1/length(ix),.5]);h.LooseInset=h.TightInset;
    h.Children(1).LineStyle='--';h.YLim(2)=0.1;h.XLim(1)=0;h.XTick=0:.2:1;
    legend(h,'$$\Gamma_N$$','$$\Gamma_G$$','Interpreter','latex','Box','off','Location','best');
    annotation('textbox',[h.Position(1),h.Position(2),0.1,0.1],'String',['$$x=',num2str(x_s(i)),'$$'],'FitBoxToText','on','Interpreter','latex','EdgeColor','none','FontSize',12);
    annotation('textbox',[h.Position(1),min(h.Position(2)+h.Position(4)+0.02,1),0,0],'String',['(',char(96+i),')'],'FitBoxToText','on','EdgeColor','none','FontSize',12);
end
temp_var(:,:,1)=n;temp_var(:,:,2)=g;temp_title={'N','G'};xlabel_pos=[0.0025,0.035];
for i=1:2
    h=draw_plot({h_s(i+2),(1:nts)*tau*ntp,temp_var(:,:,i)},'','t',temp_title{i},'Box','on','FontSize',12,'XTick',0:.2:1);%,'OuterPosition',[0,0,.5,.5]);h.LooseInset=h.TightInset;
    % h.Title.Units='normalized';h.Title.Position(2)=h.Title.Position(2)-0.02;
    h.Children(1).Color=h.ColorOrder(4,:);h.XLim(1)=0;
    annotation('textbox',[h.Position(1),h.Position(2)+h.Position(4)+0.02,0,0],'String',['(',char(98+i),')'],'FitBoxToText','on','EdgeColor','none','FontSize',12);
    h.YLabel.Position(1)=h.YLabel.Position(1)+0.04;h.XLabel.Position(2)=h.XLabel.Position(2)+xlabel_pos(i);
    % legend(cellstr(num2str(ix','$$x=%g$$')),'Interpreter','latex','Box','off','Location','best');h.XLabel.String='t';h.XLabel.Interpreter='latex';
end
legend(h,cellstr(num2str(ix','$$x=%g$$')),'Interpreter','latex','Box','off','Location','best');%h.XLabel.String='t';h.XLabel.Interpreter='latex';
%%
print('paper/fig_t','-dpng')
print('paper/fig_t','-depsc')
%% 图5
%% 频谱非线性时间平均
ix1=find(x>0.9,1);%最靠近x=1的位置
ix2=find(x>1.3,1);
iz=2;
nt=84;
data=load(sprintf('data/dat%4.4d.mat',nt),'phi','pei','deni');
phim1=amplitude(data.phi(ix1,2:end-1,iz));
phim2=amplitude(data.phi(ix2,2:end-1,iz));
denm1=amplitude(data.deni(ix1,2:end-1,iz));
pem1=amplitude(data.pei(ix1,2:end-1,iz));
denm2=amplitude(data.deni(ix2,2:end-1,iz));
pem2=amplitude(data.pei(ix2,2:end-1,iz));

denm1(1)=inf;
denm2(1)=inf;
close
h=subplot(1,2,1);
draw_plot({h,0:59,[phim1(1:60);denm1(1:60);phim2(1:60);denm2(1:60)]},'','m','','YScale','log','Box','on','OuterPosition',[0,0,.5,1],'FontSize',12);h.LooseInset=h.TightInset;
%%
phim1=zeros(1,floor(ny0));%phim1=zeros(1,floor(ny/2));
pem1=phim1;denm1=phim1;
phim2=phim1;
pem2=phim1;denm2=phim1;
for nt=150:nts
    data=load(sprintf('data/dat%4.4d.mat',nt));
    phim1=amplitude(xcorr(data.phi(ix1,2:end-1,iz)))+phim1;
    phim2=amplitude(xcorr(data.phi(ix2,2:end-1,iz)))+phim2;
    if (isdeltaf)
        pem=abs(fft(data.pei(ix,:,iz)+pe0(ix,:,iz)))/ny+pem;
        denm1=abs(fft(data.deni(ix,:,iz)+den0(ix,:,iz)))/ny+denm1;
    else
        pem1=amplitude(xcorr(data.pei(ix1,2:end-1,iz)),xcorr(data.pei(ix1,2:end-1,iz)))+pem1;
        denm1=amplitude(xcorr(data.deni(ix1,2:end-1,iz)))+denm1;
        pem2=amplitude(xcorr(data.pei(ix2,2:end-1,iz)))+pem2;
        denm2=amplitude(xcorr(data.deni(ix2,2:end-1,iz)))+denm2;
%         pem1=amplitude(data.pei(ix1,2:end-1,iz))+pem1;
%         denm1=amplitude(data.deni(ix1,2:end-1,iz))+denm1;
%         pem2=amplitude(data.pei(ix2,2:end-1,iz))+pem2;
%         denm2=amplitude(data.deni(ix2,2:end-1,iz))+denm2;
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
h=subplot(1,2,2);
draw_plot({h,1:20,[phim1(2:21);denm1(2:21);phim2(2:21);denm2(2:21)]},'','m','','XScale','log','YScale','log','Box','on','OuterPosition',[.5,0,.5,1],'FontSize',12);h.LooseInset=h.TightInset;
h.Children(1).LineStyle='--';h.Children(3).LineStyle='--';h.Children(1).Color=h.ColorOrder(2,:);h.Children(2).Color=h.ColorOrder(2,:);h.Children(3).Color=h.ColorOrder(1,:);
lgd=legend('$$\Phi(x=0.9)$$','$$N(x=0.9)$$','$$\Phi(x=1.3)$$','$$N(x=1.3)$$');lgd.Interpreter='latex';lgd.Box='off';
%% 改进
x_s=[0.9,1.3];
nt_start=150;
ix=arrayfun(@(x_s)find(x>x_s,1),x_s);
iz=2;
m_s=[0,1,3,10];
Nm=zeros(nts,floor(ny0/2)+1);Gm=Nm;Phim=Nm;
for nt=1:nts
    data=load(sprintf('data/dat%4.4d.mat',nt),'phi','deni','pei');
    %     flux_pt(:,nt+1-nt_start)=data.pei(ix,iy,iz).*data.vex(ix,iy,iz);
    %     flux_nt(:,nt+1-nt_start)=data.deni(ix,iy,iz).*data.vex(ix,iy,iz);
    Nm(nt,:)=amplitude(data.deni(ix(1),2:end-1,iz));
    Gm(nt,:)=amplitude(data.pei(ix(1),2:end-1,iz));
    Phim(nt,:)=amplitude(data.phi(ix(1),2:end-1,iz));
end
h=subplot(1,2,1);
draw_plot({h,(1:nt)*tau*ntp,Phim(:,m_s+1)},'','t','|\Phi_m|','YScale','log','Box','on','OuterPosition',[0,0,.5,1],'FontSize',12,'XTick',0:.2:1);h.LooseInset=h.TightInset;
h.Children(1).Color=h.ColorOrder(4,:);h.XLim(1)=0;
lgd=legend(cellstr(num2str(m_s')),'Box','off','Location','best','Interpreter','latex');lgd.Title.String='$m$';
annotation('textbox',[h.Position(1),h.Position(4),0.1,0.1],'String',num2str(x_s(1),'$$x=%g$$'),'FitBoxToText','on','Interpreter','latex','EdgeColor','none','FontSize',12);
%% 左图用RMS代替
plot_variables={'dPhi','dG','dN','dvx'}';%遍历顺序为先列后行
plot_temp={'phi','pei','deni','vex'};
plot_titles={'$\textrm{RMS }\delta\Phi$','$\textrm{RMS }\delta G$','$\textrm{RMS }\delta N$','$\textrm{RMS }\delta v_x$'};
iz=2;
data2=struct;
% for nt=nt_start:nts
%     data=load(sprintf('data/dat%4.4d.mat',nt),'phi','pei','vex','deni','wi');
%     for i=1:numel(plot_temp)
%         data2.(plot_temp{i})(:,:,:,nt-nt_start+1)=data.(plot_temp{i});
%     end
% end
% for i=1:numel(plot_temp)%做时间平均
%     data2.(plot_temp{i})=mean(data2.(plot_temp{i}),4);
% end
% for i=1:numel(plot_variables)
%    data2.(plot_variables{i})=std(delt(data2.(plot_temp{i})(:,:,iz)),0,2); 
% end
for nt=170:170
    data=load(sprintf('data/dat%4.4d.mat',nt),'phi','pei','vex','deni','wi');
    for i=1:numel(plot_temp)
%         data2.(plot_temp{i})(:,:,:,nt-nt_start+1)=data.(plot_temp{i});
        data2.(plot_variables{i})(:,nt-nt_start+1)=std(delt(data.(plot_temp{i})(:,:,iz)),1,2); 
    end
    data2.dn(:,nt-nt_start+1)=std(delt(data.deni(:,:,iz).*repmat(x',[1,size(data.deni,2)]).^4),1,2);
    data2.N(:,nt-nt_start+1)=mean(data.deni(:,:,iz),2);
    data2.Phi(:,nt-nt_start+1)=mean(data.phi(:,:,iz),2);%平均值
    data2.dfluxN(:,nt-nt_start+1)=std(delt(data.deni(:,:,iz).*data.vex(:,:,iz)),1,2);
end
for i=1:numel(plot_temp)%做时间平均
    data2.(plot_variables{i})=mean(data2.(plot_variables{i}),2);
end
data2.dn=mean(data2.dn,2);
data2.dfluxN=mean(data2.dfluxN,2);
h=subplot(2,2,1);
% draw_plot({h,x,cell2mat(struct2cell(rmfield(data2,'dPhi'))')},'','x','','Box','on','FontSize',12,'OuterPosition',[0,1/2,1/2,1/2]);h.LooseInset=h.TightInset;
% draw_plot({h,x,data2.dN},'$\textrm{RMS }\delta N$','x','','Box','on','FontSize',12,'OuterPosition',[0,1/2,1/2,1/2]);h.LooseInset=h.TightInset;
draw_plot({h,x,[data2.dN]},'','','\delta N','Box','on','FontSize',12,'OuterPosition',[0,1/2,1/2,1/2]);h.LooseInset=h.TightInset;%legend('$N$','$G$','Interpreter','latex','Box','off','Location','best');
h.Children(1).Color=h.ColorOrder(2,:);
h=subplot(2,2,3);
% draw_plot({h,x,data2.dPhi},'$\textrm{RMS }\delta\Phi$','x','','Box','on','FontSize',12,'OuterPosition',[0,0,1/2,1/2]);h.LooseInset=h.TightInset;
draw_plot({h,x,[data2.dPhi]},'','x','\delta\Phi','Box','on','FontSize',12,'OuterPosition',[0,0,1/2,1/2]);h.LooseInset=h.TightInset;%legend('$\Phi$','$\Gamma_N$','Interpreter','latex','Box','off','Location','best');

%% 右图
phim=zeros(length(x_s),floor(ny0/2)+1);pem=phim;denm=phim;
for nt=nt_start:nts
    data=load(sprintf('data/dat%4.4d.mat',nt),'phi','deni','pei');
    for i=1:length(x_s)
        phim(i,:)=amplitude(data.phi(ix(i),2:end-1,iz))+phim(i,:);
        denm(i,:)=amplitude(data.deni(ix(i),2:end-1,iz))+denm(i,:);
        pem(i,:)=amplitude(data.pei(ix(i),2:end-1,iz))+pem(i,:);
    end
end
phim=phim/(nts-nt_start+1);pem=pem/(nts-nt_start+1);denm=denm/(nts-nt_start+1);
% denm(:,1)=inf;
% h=subplot(1,2,2,'NextPlot','replacechildren','LineStyleOrder',{'-','--'});h.ColorOrder(length(x_s)+1:end,:)=[];
% h=draw_plot({h,1:30,[phim(:,2:31);denm(:,2:31)]},'','m','','YScale','log','Box','on','FontSize',12,'OuterPosition',[0.5,0,0.5,1]);h.LooseInset=h.TightInset;
h=subplot(1,2,2,'NextPlot','replacechildren','LineStyleOrder',{'-','--'});h.ColorOrder(length(x_s)+1:end,:)=[];
h=draw_plot({h,1:20,[phim(:,2:21).^2;denm(:,2:21).^2]},'','m','|\delta\Phi|^2','XScale','log','YScale','log','Box','on','FontSize',12,'OuterPosition',[1/2,0,1/2,1]);h.LooseInset=h.TightInset;
% h.Children(1).LineStyle='--';h.Children(2).LineStyle='--';h.Children(1).Color=h.ColorOrder(2,:);h.Children(3).Color=h.ColorOrder(2,:);h.Children(2).Color=h.ColorOrder(1,:);
lgd=legend([cellstr(num2str(x_s','$$\\Phi(x=%g)$$'));cellstr(num2str(x_s','$$N(x=%g)$$'))],'Interpreter','latex','Box','off','Location','best');
%% 右图固定位置, 加vx
x_s=1.3;
nt_start=150;
ix=arrayfun(@(x_s)find(x>x_s,1),x_s);
phim=zeros(length(x_s),floor(ny0/2)+1);pem=phim;denm=phim;vxm=phim;fluxnm=phim;
for nt=nt_start:nts
    data=load(sprintf('data/dat%4.4d.mat',nt),'phi','deni','pei','vex');
    for i=1:length(x_s)
        phim(i,:)=amplitude(data.phi(ix(i),2:end-1,iz))+phim(i,:);
        denm(i,:)=amplitude(data.deni(ix(i),2:end-1,iz))+denm(i,:);
        pem(i,:)=amplitude(data.pei(ix(i),2:end-1,iz))+pem(i,:);
        vxm(i,:)=amplitude(data.vex(ix(i),2:end-1,iz))+vxm(i,:);
        fluxnm(i,:)=amplitude(data.vex(ix(i),2:end-1,iz).*data.deni(ix(i),2:end-1,iz))+fluxnm(i,:);
    end
end
phim=phim/(nts-nt_start+1);pem=pem/(nts-nt_start+1);denm=denm/(nts-nt_start+1);vxm=vxm/(nts-nt_start+1);fluxnm=fluxnm/(nts-nt_start+1);
h=subplot(2,2,[2,4]);
% h=draw_plot({h,1:20,[phim(:,2:21).^2;denm(:,2:21).^2;vxm(:,2:21).^2]},'','m','','XScale','log','YScale','log','Box','on','FontSize',12,'OuterPosition',[1/2,0,1/2,1]);h.LooseInset=h.TightInset;
h=draw_plot({h,1:20,[phim(:,2:21).^2;denm(:,2:21).^2;vxm(:,2:21).^2]},'','m','','XScale','log','YScale','log','Box','on','FontSize',12,'OuterPosition',[1/2,0,1/2,1]);h.LooseInset=h.TightInset;
h.Children(1).Color=h.ColorOrder(4,:);
legend({'$|\delta\Phi|^2$','$|\delta N|^2$','$|\delta v_x|^2$'},'Interpreter','latex','Box','off','Location','best');
%%
print('paper/fig_spectrum','-dpng');
print('F:/科研/paper/extended.190627/fig_spectrum','-depsc','-tiff');
%% 图6
%% eta和omegapd特定x位置随时间演化
ix=[.7,.9];
nt_start=50;
ix_index=arrayfun(@(ix)find(x>ix,1),ix);
baseline={5/3,2/3};
baseline_text={'$1.67$','$0.67$'};
plot_variables={'omegap','eta'}';%variables name, 遍历顺序为先列后行
plot_titles={'\omega_p^*/\omega_d','\eta'};
subplot_grid=num2cell(size(plot_variables'));
data2=struct;
for nt=nt_start:nts
    data=load(sprintf('data/dat%4.4d.mat',nt));
    if isdeltaf
        pe0=pe0(:,:,2)';%%pe0为平均量
        den0=den0(:,:,2)';
    else
        pe0=mean(data.pei(:,2:end-1,2),2)';
        den0=mean(data.deni(:,2:end-1,2),2)';
    end
    dpe=gradient(pe0,dx);%dG/dx
    dden=gradient(den0,dx);%dN/dx
    eta=(dpe./pe0+4*gamma./x)./(dden./den0+4./x)-1;
    data2.eta(:,nt+1-nt_start)=eta(ix_index);
    omegap=dpe./pe0/4.*x+gamma;
    data2.omegap(:,nt+1-nt_start)=omegap(ix_index);
end
close
for i=1:numel(plot_variables)
    subplot(subplot_grid{:},i);
    h=draw_plot({(nt_start:nts)*tau*ntp,data2.(plot_variables{i})},'','t',plot_titles{i},'Box','on','FontSize',12,...
        'OuterPosition',[[0,1]+[rem(i-1,size(plot_variables,1)),-1-floor((i-1)/size(plot_variables,1))]./size(plot_variables),1./size(plot_variables)]);
%     h.YLabel.Position(1)=h.YLabel.Position(1)+0.07;
    h.LooseInset=h.TightInset;
    h.Children(1).LineStyle='--';
    %     line(h,'XData',h.XLim,'YData',baseline{i}*ones(1,2),'LineStyle','-.')
    %     [h.YTick,sortIndex]=sort([baseline{i},h.YTick]);h.YTickLabel{find(sortIndex==1,1)}=baseline_text{i};h.TickLabelInterpreter='latex';
    legend(cellstr(num2str(ix','$$x=%g$$')),'Interpreter','latex','Box','off','Location','best');
end
print('paper/fig_eta','-dpng');
print('paper/fig_eta','-deps');


%% 图5
%% 环向平均N随时间演化, 给定x
ix=[0.8,0.9,1];
ix_index=arrayfun(@(ix)find(x>ix,1),ix);
n=zeros(length(ix),nts);g=n;
for nt=1:nts
    data=load(sprintf('data/dat%4.4d.mat',nt),'deni','pei');
    n(:,nt)=mean(data.deni(ix_index,2:end-1,2),2);
    g(:,nt)=mean(data.pei(ix_index,2:end-1,2),2);
end
%%
close;
h=subplot(1,2,1);
draw_plot({(1:nts)*tau*ntp,n},'$$N$$','t','','Box','on','FontSize',12,'OuterPosition',[0,0,.5,1]);h.LooseInset=h.TightInset;
% h.Title.Units='normalized';h.Title.Position(2)=h.Title.Position(2)-0.02;
h.Children(1).Color='#7E2F8E';
legend(cellstr(num2str(ix','$$x=%g$$')),'Interpreter','latex','Box','off','Location','best');
h=subplot(1,2,2);
% draw_plot({(1:nts)*tau*ntp,(n).*repmat(ix',[1,nts]).^4},'$$n$$','t','','Box','on','FontSize',12,'OuterPosition',[(2-1)/2,0,1/2,1]);h.LooseInset=h.TightInset;
draw_plot({(1:nts)*tau*ntp,g},'$$G$$','t','','Box','on','FontSize',12,'OuterPosition',[0.5,0,0.5,1]);h.LooseInset=h.TightInset;
% h.Title.Units='normalized';h.Title.Position(2)=h.Title.Position(2)-0.02;
h.Children(1).Color='#7E2F8E';
legend(cellstr(num2str(ix','$$x=%g$$')),'Interpreter','latex','Box','off','Location','best');
% h=subplot(2,2,3);
% draw_plot({(1:nts)*tau*ntp,(n).*repmat(ix',[1,nts]).^4},'$$n$$','t','','Box','on','FontSize',12,'OuterPosition',[0,0,.5,.5]);h.LooseInset=h.TightInset;
% % h.Title.Units='normalized';h.Title.Position(2)=h.Title.Position(2)-0.02;
% h.Children(1).Color='#7E2F8E';
% legend(cellstr(num2str(ix','$$x=%g$$')),'Interpreter','latex','Box','off','Location','best');
% h=subplot(2,2,4);
% draw_plot({(1:nts)*tau*ntp,(g).*repmat(ix',[1,nts]).^(4*gamma)},'$$p$$','t','','Box','on','FontSize',12,'OuterPosition',[0.5,0,.5,.5]);h.LooseInset=h.TightInset;
% % draw_plot({(1:nts)*tau*ntp,g},'$$G$$','t','','Box','on','FontSize',12,'OuterPosition',[(2-1)/2,0,1/2,1]);h.LooseInset=h.TightInset;
% % h.Title.Units='normalized';h.Title.Position(2)=h.Title.Position(2)-0.02;
% h.Children(1).Color='#7E2F8E';
% legend(cellstr(num2str(ix','$$x=%g$$')),'Interpreter','latex','Box','off','Location','best');
print('paper/fig_N_t','-dpng')
print('paper/fig_N_t','-depsc','-tiff')
%% 图4
%% 随时间演化的flux
x_s=[0.8,0.9];
ix=arrayfun(@(x_s)find(x>x_s,1),x_s);
iy=ny/2;
iz=2;
nt_start=1;
flux_pt=zeros(length(ix),nts-nt_start+1);
flux_nt=flux_pt;
for nt=nt_start:nts
    data=load(sprintf('data/dat%4.4d.mat',nt),'pei','vex','deni');
    %     flux_pt(:,nt+1-nt_start)=data.pei(ix,iy,iz).*data.vex(ix,iy,iz);
    %     flux_nt(:,nt+1-nt_start)=data.deni(ix,iy,iz).*data.vex(ix,iy,iz);
    flux_pt(:,nt+1-nt_start)=mean(data.pei(ix,2:end-1,iz).*data.vex(ix,2:end-1,iz),2);
    flux_nt(:,nt+1-nt_start)=mean(data.deni(ix,2:end-1,iz).*data.vex(ix,2:end-1,iz),2);
end
for nt=nt_start:nts
    data=load(sprintf('data/dat%4.4d.mat',nt),'pei','vex','deni');
    %     flux_pt(:,nt+1-nt_start)=data.pei(ix,iy,iz).*data.vex(ix,iy,iz);
    %     flux_nt(:,nt+1-nt_start)=data.deni(ix,iy,iz).*data.vex(ix,iy,iz);
    flux_pt(:,nt+1-nt_start)=data.pei(ix,iy,iz).*data.vex(ix,iy,iz);
    flux_nt(:,nt+1-nt_start)=data.deni(ix,iy,iz).*data.vex(ix,iy,iz);
end
%%
close
for i=1:length(ix)
    h=subplot(1,length(ix),i);
    draw_plot({h,(nt_start:nts)*tau*ntp,[flux_pt(i,:);flux_nt(i,:)]},'','t','','Box','on','FontSize',12,'OuterPosition',[(i-1)/length(ix),0,1/length(ix),1]);
    h.LooseInset=h.TightInset;
    h.Children(2).LineStyle='--';
    legend('$$\Gamma_G$$','$$\Gamma_N$$','Interpreter','latex','Box','off','Location','best');
    annotation('textbox',[h.Position(1),h.Position(4),0.1,0.1],'String',['$$x=',num2str(x_s(i)),'$$'],'FitBoxToText','on','Interpreter','latex','EdgeColor','none','FontSize',12);
end
print('paper/fig_flux_t_2','-dpng');
print('paper/fig_flux_t_2','-depsc','-tiff');
function vm=amplitude(v)
% AMPLITUDE calculate the amplitude of phase
% v is the real vector
% vm has length of floor(n/2)+1, with phase number 0:floor(n/2)+1
n=length(v);
vm=abs(fft(v))/n;
vm(2:floor((n+1)/2))=2*vm(2:floor((n+1)/2));
vm(floor(n/2)+2:end)=[];
end