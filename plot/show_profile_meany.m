% 环向平均profile
close all
nt_slice=[1,16,30];%初始状态, 线性阶段, 非线性阶段
subplot_grid={3,2};
if (isdeltaf)
plot_variables={'pei(:,2:ny0+1,iz)+pe0','p',...%variables name
    'deni(:,2:ny0+1,iz)+den0','n',...
    'wi','phi'};
else
    plot_variables={'pei','p',...%variables name
    'deni','n',...
    'wi','phi'};
end
plot_titles={'G','p_e','N','\langle n\rangle','w','\Phi'};%variables title
plot_subtitles={'in initial state','in linear state','in nonlinear state'};
for i=1:length(nt_slice)
    clf;
    load(sprintf('data/dat%4.4d.mat',i))
    if (isdeltaf)
        n=(deni+den0).*repmat(x',[1,size(deni,2),size(deni,3)]).^4;
        p=(pei+pe0).*repmat(x',[1,size(deni,2),size(deni,3)]).^(4*gamma);
    else
        n=deni.*repmat(x',[1,size(deni,2),size(deni,3)]).^4;
        p=pei.*repmat(x',[1,size(deni,2),size(deni,3)]).^(4*gamma);
    end
    for j=1:length(plot_variables)
        subplot(subplot_grid{:},j);
        draw_plot({x,rms(eval([plot_variables{j},'(:,2:ny0+1,iz)']),2)},['$$',plot_titles{j},'$$'],'x','');
    end
    suptitle({'profile with $$y$$ rms average',plot_subtitles{i}});
    print(['plot/mean_y_',num2str(i)],'-dpng');
end
close