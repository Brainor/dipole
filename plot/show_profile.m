%% 环向和时间平均profile
nt_s=150:200;
plot_variables={'pe0','p';...%variables name
    'den0','n';...
    'phi','wi'}';%遍历顺序为先列后行
plot_titles={'G','p_e','N','\langle n\rangle','\Phi','w'};
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
    data.n=mean(data.den0.*repmat(x',[1,size(deni,2)]).^4,2);
    data.p=mean(data.pe0.*repmat(x',[1,size(deni,2)]).^(4*gamma),2);
    for i=1:numel(plot_variables)
        data2.(plot_variables{i})(:,nt-nt_s(1)+1)=data.(plot_variables{i});
    end
end
for i=1:numel(plot_variables)%做时间平均
    data2.(plot_variables{i})=mean(data2.(plot_variables{i}),2);
end
if(sp); data2.pe0=[data2.pe0';s_p/tau];end
if(sn); data2.den0=[data2.den0';s_den/tau];end%加上源
% 横坐标为x
for i=1:numel(plot_variables)
    subplot(subplot_grid{:},i);
    draw_plot({x,data2.(plot_variables{i})},['$$',plot_titles{i},'$$'],'x','','Box','on');
end
suptitle({'profile with time and $$y$$ average in nonlinear phase',['slice ',sprintf('%d-%d',nt_s(1),nt_s(end))]})
print('plot/profile_x_meanyt','-dpng');
clf;
% 横坐标为L
for i=1:numel(plot_variables)
    subplot(subplot_grid{:},i);
    draw_plot({1./x,data2.(plot_variables{i})},['$$',plot_titles{i},'$$'],'x','');
end
suptitle({'profile with time and $$y$$ average in nonlinear phase','in radial direction'})
print('plot/profile_L_meanyt','-dpng');
clf;
close

%% 环向平均profile
nt_slice=[1,72,157];%初始状态, 线性阶段, 非线性阶段
plot_variables={'pe0','den0';
    'fluxp','fluxn';
    'p','n'}';%variables name, 遍历顺序为先列后行
plot_titles={'\textrm{RMS }\delta\Phi','\textrm{RMS }\delta G','G_0','N_0','\Gamma_G','\Gamma_N','p_e','\langle n\rangle'};
plot_variables={'pe0','den0';
    'fluxp','fluxn';
    'p','n'}';%variables name, 遍历顺序为先列后行
plot_titles={'G','N','\Gamma_G','\Gamma_N','p_e','n'};
plot_subtitles={'in initial state','in linear state','in nonlinear state'};
subplot_grid=num2cell(size(plot_variables'));
figure
for i=1:length(nt_slice)
    clf;
    data=load(sprintf('data/dat%4.4d.mat',nt_slice(i)),'deni','phi','vex','pei');
    if (isdeltaf)
        data.pe0=mean(pe0(:,:,iz)-bc_p,2);
        data.den0=mean(den0(:,:,iz)-bc_n,2);
    else
        data.pe0=mean(data.pei(:,:,iz)-bc_p,2);
        data.den0=mean(data.deni(:,:,iz)-bc_n,2);
    end
    data.phi=std(delt(data.phi(:,:,iz)),0,2);
    data.dpe=std(delt(data.pei(:,:,iz)),0,2);
    data.fluxp=mean(data.vex(:,:,iz).*data.pei(:,:,iz),2);
    data.fluxn=mean(data.vex(:,:,iz).*data.deni(:,:,iz),2);
    data.n=mean(data.den0.*repmat(x',[1,size(deni,2)]).^4,2);
    data.p=mean(data.pe0.*repmat(x',[1,size(deni,2)]).^(4*gamma),2);
    if(sp); data.pe0=[data.pe0';s_p/tau];end
    if(sn); data.den0=[data.den0';s_den/tau];end%加上源
    for j=1:numel(plot_variables)
        subplot(subplot_grid{:},j);
        draw_plot({x,data.(plot_variables{j})},['$$',plot_titles{j},'$$'],'x','','Box','on');
    end
    suptitle({'profile with $$y$$ mean average',plot_subtitles{i},num2str(nt_slice(i))});
    print(['plot/profile_x_meany_',num2str(i)],'-dpng');
end
close
%% flux随时间演化, movie
plot_variables={'pe0','den0';
    'fluxp','fluxn';
    'p','n'}';%variables name, 遍历顺序为先列后行
plot_titles={'G','N','Gv_x','Nv_x','p_e','\langle n\rangle'};
subplot_grid=num2cell(size(plot_variables'));
v = VideoWriter('flux','MPEG-4');
v.FrameRate=10;
tic
mkdir('temp')
for nt=1:nts
    data=load(sprintf('data/dat%4.4d.mat',nt),'pei','vex','deni');
    if (isdeltaf)
        data.pe0=mean(pe0(:,:,iz)-bc_p,2);
        data.den0=mean(den0(:,:,iz)-bc_n,2);
    else
        data.pe0=mean(data.pei(:,:,iz)-bc_p,2);
        data.den0=mean(data.deni(:,:,iz)-bc_n,2);
    end
    data.fluxp=mean(data.vex(:,:,iz).*data.pei(:,:,iz),2);
    data.fluxn=mean(data.vex(:,:,iz).*data.deni(:,:,iz),2);
    data.n=mean(data.den0.*repmat(x',[1,size(deni,2)]).^4,2);
    data.p=mean(data.pe0.*repmat(x',[1,size(deni,2)]).^(4*gamma),2);
    if(sp); data.pe0=[data.pe0';s_p/tau];end
    if(sn); data.den0=[data.den0';s_den/tau];end%加上源
    for i=1:numel(plot_variables)
        subplot(subplot_grid{:},i);
        draw_plot({x,data.(plot_variables{i})},['$$',plot_titles{i},'$$'],'x','');
    end
    suptitle(num2str(nt),'FontSize',8);
    print(sprintf('temp/%4.4d',nt),'-dpng');
end
makevideo('temp/*.png','flux');
rmdir temp s
toc
%% 环向平均N随时间演化, 给定x
ix=[0.7,0.8,0.9,1];
ix_index=arrayfun(@(ix)find(x>ix,1),ix);
n=zeros(length(ix),nts);
for nt=1:nts
    data=load(sprintf('data/dat%4.4d.mat',nt),'deni');
    n(:,nt)=mean(data.deni(ix_index,:,2),2);
end
subplot(1,2,1)
draw_plot({1:nts,n},'density profile $$N$$ with toroidal average','t','','Box','off');
subplot(1,2,2)
draw_plot({1:nts,(n-1).*repmat(ix',[1,nts]).^4},'density profile $$n$$ with toroidal average','t','','Box','off');
legend(cellstr(num2str(ix','$$x=%g$$')),'Interpreter','latex','Box','off');
print('plot/N_t','-dpng')
%% 环向平均flux随时间演化, 给定x
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
subplot(1,2,1)
draw_plot({(nt_start:nts),[flux_pt(1,:);flux_nt(1,:)]},'flux time evolution at $$x=0.9$$','t','','Box','on');
legend('$$\Gamma_G$$','$$\Gamma_N$$','Interpreter','latex','Box','off');
subplot(1,2,2)
draw_plot({(nt_start:nts),[flux_pt(2,:);flux_nt(2,:)]},'flux time evolution at $$x=0.7$$','t','','Box','on');
legend('$$\Gamma_G$$','$$\Gamma_N$$','Interpreter','latex','Box','off');
suptitle('slice 1-200');
print('plot/flux_ts','-dpng');
%% Phi和N的关联函数
ix=find(x>0.9,1);
iy=[2,ny/2];%ix和iy只能有一个量为向量, 需要改subplot(1,3,3)的title
iz=2;
nt_s=1:200;
plot_variables={'phi','pei','deni','vex',}';
plot_titles={'\Phi','G','N','v_x'};%variables title
data2=struct;
for nt=nt_s
    data=load(sprintf('data/dat%4.4d.mat',nt),'phi','pei','vex','deni');
    for i=1:numel(plot_variables)
        data2.(plot_variables{i})(:,nt-nt_s(1)+1)=data.(plot_variables{i})(ix,iy,iz);
    end
end
for i=1:numel(plot_variables)%减去时间平均项
    data2.(plot_variables{i})=delt(data2.(plot_variables{i}),2);
end
%互相关
cross=zeros(numel(plot_variables)-1,2*length(nt_s)-1);
for i=2:numel(plot_variables)
    [cross(i-1,:),lags]=xcorr(data2.phi(1,:),data2.(plot_variables{i})(1,:));
end
lags=lags*tau*ntp;
subplot(1,3,1)
draw_plot({lags,cross},'Cross-Correlation of variables and $$\Phi$$','t');
legend(cellfun(@(x)['$$',x,'$$'],plot_titles(2:end),'UniformOutput',false),'Interpreter','latex','Location','best','Box','off');
%自相关
cross=zeros(numel(plot_variables),2*length(nt_s)-1);
for i=1:numel(plot_variables)
    cross(i,:)=xcorr(data2.(plot_variables{i})(1,:));
end
subplot(1,2,1)
draw_plot({lags,cross},'Auto-Correlation of variables','t','','Box','on');
legend(cellfun(@(x)['$$',x,'$$'],plot_titles,'UniformOutput',false),'Interpreter','latex','Location','best','Box','off');
%相同参数不同位置互相关
cross=zeros(numel(plot_variables),2*length(nt_s)-1);
for i=1:numel(plot_variables)
    cross(i,:)=xcorr(data2.(plot_variables{i})(1,:),data2.(plot_variables{i})(2,:));
end
subplot(1,2,2)
draw_plot({lags,cross},{'Cross-Correlation of variables','at $$y=0$$ and $$y=\pi$$'},'t','','Box','on');
legend(cellfun(@(x)['$$',x,'$$'],plot_titles,'UniformOutput',false),'Interpreter','latex','Location','best','Box','off');

print('plot/correlation','-dpng');