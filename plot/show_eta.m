%% 做eta和omega_p/omega_d
close all
plot_variables={'eta','omegap'}';%variables name, 遍历顺序为先列后行
plot_titles={'\eta','\frac{\omega_p}{\omega_d^\star}'};
subplot_grid=size(plot_variables');
subplot_gridcell=num2cell(subplot_grid);
v = VideoWriter('eta','MPEG-4');
v.FrameRate=10;
open(v);
for nt=1:nts
    load(sprintf('data/dat%4.4d.mat',nt))
    if isdeltaf
        pe0=pe0(:,:,2)';%%pe0为平均量
        den0=den0(:,:,2)';
    else
        pe0=mean(pei(:,:,2),2)';
        den0=mean(deni(:,:,2),2)';
    end
    dpe=gradient(pe0,dx);%dG/dx
    dden=gradient(den0,dx);%dN/dx
    eta=(dpe./pe0+4*gamma./x)./(dden./den0+4./x)-1;
    omegap=dpe./pe0/4.*x+gamma;
    for i=1:length(plot_variables)
        subplot(subplot_gridcell{:},i);
        draw_plot({x,eval(plot_variables{i})},['$$',plot_titles{i},'$$'],'x','');
    end
    suptitle(num2str(nt),'FontSize',8);
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);