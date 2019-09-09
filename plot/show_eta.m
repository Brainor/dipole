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
    for i=1:numel(plot_variables)
        subplot(subplot_gridcell{:},i);
        draw_plot({x,eval(plot_variables{i})},['$$',plot_titles{i},'$$'],'x','');
    end
    suptitle(num2str(nt),'FontSize',8);
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);
%% eta和omegapd特定x位置随时间演化
ix=[0.7,0.9];
ix_index=arrayfun(@(ix)find(x>ix,1),ix);
baseline={5/3,2/3};
baseline_text={'1.67','0.67'};
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
clf
for i=1:numel(plot_variables)
    subplot(subplot_grid{:},i);
    h=draw_plot({1:nts,data2.(plot_variables{i})},['$$',plot_titles{i},'$$'],'t','','Box','on');
    line(h,'XData',h.XLim,'YData',baseline{i}*ones(1,2),'LineStyle','--')
    h.YTick=sort([baseline{i},h.YTick]);h.YTickLabel{2}=baseline_text{i};
    legend(cellstr(num2str(ix','$$x=%g$$')),'Interpreter','latex','Box','off');
end
print('plot/eta','-dpng');