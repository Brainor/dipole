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
draw_plot({(nt_start:nts),[flux_pt;flux_nt]},['flux time evolution at slice ',num2str(nt_start),'-',num2str(nts)],'t','','Box','on');
legend('$$\Gamma_G(x=0.9)$$','$$\Gamma_N(x=0.9)$$','$$\Gamma_G(x=0.5)$$','$$\Gamma_N(x=0.5)$$','Interpreter','latex','Box','off');
print('plot/flux_t','-dpng');
%% temp
subplot(1,2,1)
draw_plot({(nt_start:nts),[flux_pt(1,:);flux_nt(1,:)]},'flux time evolution at slice 1-200','t','','Box','on');
legend('$$\Gamma_G(x=0.9)$$','$$\Gamma_N(x=0.9)$$','Interpreter','latex','Box','off');
subplot(1,2,2)
draw_plot({(nt_start:nts),[flux_pt(2,:);flux_nt(2,:)]},'flux time evolution at slice 1-200','t','','Box','on');
legend('$$\Gamma_G(x=0.7)$$','$$\Gamma_N(x=0.7)$$','Interpreter','latex','Box','off');
print('plot/flux_ts','-dpng');
%% 频数图
[counts,edges]=histcounts(flux_nt(120:end),8);
draw_plot({edges(1:end-1),counts},'PDF at 120-200');

%% 时间movie
plot_variables={'mean(data.pei(:,:,iz),2)','mean(data.deni(:,:,iz),2)';
    'mean(data.vex(:,:,iz).*data.pei(:,:,iz),2)','mean(data.vex(:,:,iz).*data.deni(:,:,iz),2)'}';%variables name, 遍历顺序为先列后行
if(sp); plot_variables{1}='[mean(data.pei(:,:,iz),2)'';s_p/tau+bc_p(1)]';end
if(sn); plot_variables{2}='[mean(data.deni(:,:,iz),2)'';s_den/tau+bc_n(1)]';end%加上源
plot_titles={'\langle G\rangle','\langle N\rangle','v_xG','v_xN','n','p_e'};
subplot_grid=size(plot_variables');
subplot_gridcell=num2cell(subplot_grid);
v = VideoWriter('flux','MPEG-4');
v.FrameRate=10;
open(v);
tic
for nt=1:nts
    data=load(sprintf('data/dat%4.4d.mat',nt),'pei','vex','deni');
    for i=1:numel(plot_variables)
        subplot(subplot_gridcell{:},i);
        draw_plot({x,eval(plot_variables{i})},['$$',plot_titles{i},'$$'],'x','');
    end
    suptitle(num2str(nt),'FontSize',8);
    frame = getframe(gcf);
    writeVideo(v,frame);
end
toc
close(v);
%% 时间movie.测试
plot_variables={'pe0','den0';
    'fluxp','fluxn'}';%variables name, 遍历顺序为先列后行
plot_titles={'\langle G\rangle','\langle N\rangle','v_xG','v_xN'};
subplot_grid=num2cell(size(plot_variables'));
v = VideoWriter('flux','MPEG-4');
v.FrameRate=10;
tic
mkdir('temp')
parfor nt=1:nts
    data=load(sprintf('data/dat%4.4d.mat',nt),'pei','vex','deni');
    data.pe0=[mean(data.pei(:,:,iz),2)';s_p/tau+bc_p(1)];
    data.den0=[mean(data.deni(:,:,iz),2)';s_den/tau+bc_n(1)];
    data.fluxp=mean(data.vex(:,:,iz).*data.pei(:,:,iz),2);
    data.fluxn=mean(data.vex(:,:,iz).*data.deni(:,:,iz),2);
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
%% 用线性方式求解flux
nt=157;
ix=find(x>1,1);%最靠近x=1的位置
data=load(sprintf('data/dat%4.4d.mat',nt),'pei','vex','deni','phi');
data.pe0=[mean(data.pei(:,2:end-1,iz),2)';s_p/tau+bc_p(1)];
data.den0=[mean(data.deni(:,2:end-1,iz),2)';s_den/tau+bc_n(1)];
data.fluxG=mean(data.vex(:,2:end-1,iz).*data.pei(:,2:end-1,iz),2);
data.fluxN=mean(data.vex(:,2:end-1,iz).*data.deni(:,2:end-1,iz),2);

phim=fft(data.phi(:,2:end-1,iz),[],2)/ny0;
pem=fft(data.pei(:,2:end-1,iz),[],2)/ny0;
denm=fft(data.deni(:,2:end-1,iz),[],2)/ny0;
for i=1:ny/2
    data.fluxN2(:,i)=2*real(1i*(i-1)*phim(:,i).*conj(denm(:,i)));
    data.fluxG2(:,i)=2*real(1i*(i-1)*phim(:,i).*conj(pem(:,i)));
    data.fluxG3(:,i)=-i*abs(phim(:,i)).*abs(pem(:,i)).*sin(angle(phim(:,9))-angle(pem(:,9)));
end
subplot(2,2,1)
draw_plot({x,[sum(data.fluxG2,2),data.fluxG]},'$$\Gamma_G$$','x');
subplot(2,2,2)
draw_plot({x,[sum(data.fluxN2,2),data.fluxN]},'$$\Gamma_N$$','x');
subplot(2,2,1)
draw_plot({x,[data.fluxG2(:,3),data.fluxG]},'$$\Gamma_G$$','x');
subplot(2,2,2)
draw_plot({x,[data.fluxN2(:,3),data.fluxN]},'$$\Gamma_N$$','x');
legend('$$\Re\{im\delta\Phi_m\delta G_m^*\}$$','$$\widetilde{G}v_x$$','Interpreter','latex','Location','best','Box','off');
print('plot\test','-dpng');

subplot(2,2,3)
draw_plot({x,[angle(phim(:,2)),angle(pem(:,2)),angle(denm(:,2))]},'phase($$m=3$$)','x','\varphi');
legend('$$\Phi$$','$$G$$','$$N$$','Interpreter','latex','Box','off');
subplot(2,2,4)
draw_plot({x,[sin(angle(phim(:,3))-angle(pem(:,3))),sin(angle(phim(:,3))-angle(denm(:,3)))]},'phase difference($$m=2$$)','x','\Delta\varphi');
legend('$$\theta_\Phi-\theta_G$$','$$\theta_\Phi-\theta_N$$','Interpreter','latex','Box','off');
print('plot\test','-dpng');
%%
data=load(sprintf('data/dat%4.4d.mat',nt),'pei','vex','deni','phi');
dpe=delt(data.pei(ix,:,iz),2);
dphi=delt(data.phi(ix,:,iz),2);
dden=delt(data.deni(ix,:,iz),2);
draw_plot({-dy:dy:aly,[dpe;dphi;dden]},'$$x=1$$','y');
legend('$$\widetilde{G}$$','$$\widetilde{\Phi}$$','$$\widetilde{N}$$','Interpreter','latex','Box','off');
print('plot\test','-dpng');