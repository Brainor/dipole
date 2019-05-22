close all
global x y
%% draw the time evolution
iz=2;%2d
y=-dy:dy:aly;
subplot_grid=[4,3];%subplot的前两个参量
if (isdeltaf)%设定作图的变量
    plot_variables={'phi','pei','deni',...%变量名称
        'wi','pei+pe0','deni+den0',...
        'vey','p','n',...
        '(pei+pe0).*vex','p.*vex','n.*vex'};
else
    %     deni=deni-1;%因为N=ndV+1
    plot_variables={'phi','sp0(pei)','sp0(deni)',...%变量名称
        'wi','pei','deni',...
        'vey','p','n',...
        'pei.*vex','p.*vex','n.*vex'};
end
plot_titles={'\Phi','\widetilde{G}','\widetilde{N}','w','G','N','v_y','p_e','n','Gv_x','p_e v_x','nv_x'};%subplot标题
for nt=1:nts
    load(sprintf('data/dat%4.4d.mat',nt))
    if (isdeltaf)%设定作图的变量
        n=(deni+den0).*repmat(x',[1,size(deni,2),size(deni,3)]).^4;
        p=(pei+pe0).*repmat(x',[1,size(deni,2),size(deni,3)]).^(4*gamma);
    else
        n=deni.*repmat(x',[1,size(deni,2),size(deni,3)]).^4;
        p=pei.*repmat(x',[1,size(deni,2),size(deni,3)]).^(4*gamma);
    end
    pvx=p.*vex;
    T=p./n;
    for i=1:length(plot_variables)
        draw_pcolor([subplot_grid,i],{eval(plot_variables{i}),y,x},['$$',plot_titles{i},'$$'],'y','x');
    end
    suptitle(num2str(nt),'FontSize',8);
    print(sprintf('plot/%4.4d',nt),'-dpng')
    clf;
end
makevideo('plot/0*.png','time_evolution');
%% profile in the last slice, w/o average
subplot_grid={4,2};
plot_variables={'pei','p',...%variables name
    'deni','n',...
    'phi','wi',...
    'vey','pvx'};
plot_titles={'G','p_e','N','\langle n\rangle','\Phi','w','v_y','p_e v_x'};%variables title
n=deni.*repmat(x',[1,size(deni,2),size(deni,3)]).^4;
p=pei.*repmat(x',[1,size(deni,2),size(deni,3)]).^(4*gamma);
for i=1:length(plot_variables)
    subplot(subplot_grid{:},i);
    draw_plot({x,eval([plot_variables{i},'(:,floor(ny/2),iz)'])},['$$',plot_titles{i},'$$'],'x','');
end
suptitle('last profile at $$y=\pi$$')
print('plot/prof_last','-dpng');
close