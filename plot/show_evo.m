close all
%% draw the time evolution
tic
iz=2;%2d
y=-dy:dy:aly;
plot_variables={'phi','pei','deni';...%变量名称
        'wi','G','N';...
        'vex','fluxG','fluxN';...
        'vey','fluxp','fluxn'}';%遍历顺序为先列后行
plot_titles={'\Phi','\widetilde{G}','\widetilde{N}','w','G','N','v_x','\widetilde{G}v_x','\widetilde{N}v_x','v_y','\widetilde{p}_e v_x','\widetilde{n}v_x'};
subplot_grid=num2cell(size(plot_variables'));
for nt=1:nts
    data=load(sprintf('data/dat%4.4d.mat',nt),'phi','vex','deni','pei','wi','vey');
    if (isdeltaf)%设定作图的变量
        data.N=data.deni+den0;
        data.G=data.pei+pe0;
    else
        data.N=data.deni;data.deni=delt(data.deni);
        data.G=data.pei;data.pei=delt(data.pei);
    end
    data.p=data.G.*repmat(x',[1,size(deni,2),size(deni,3)]).^(4*gamma);
    data.n=data.N.*repmat(x',[1,size(deni,2),size(deni,3)]).^4;
    data.fluxG=data.vex.*data.pei;
    data.fluxN=data.vex.*data.deni;
    data.fluxp=data.vex.*delt(data.p,2);
    data.fluxn=data.vex.*delt(data.n,2);
    for i=1:numel(plot_variables)
        draw_pcolor([subplot_grid(:)',{i}],{data.(plot_variables{i}),y,x},['$$',plot_titles{i},'$$'],'y','x');
    end
    suptitle(num2str(nt),'FontSize',8);
    print(sprintf('plot/%4.4d',nt),'-dpng')
    clf;
end
toc
makevideo('plot/0*.png','time_evolution');
%% profile in the last slice, w/o average
subplot_grid={4,2};
plot_variables={'pei','p',...%variables name
    'deni','n',...
    'phi','wi',...
    'vey','fluxp'};
plot_titles={'G','p_e','N','\langle n\rangle','\Phi','w','v_y','p_e v_x'};%variables title
for i=1:length(plot_variables)
    subplot(subplot_grid{:},i);
    draw_plot({x,data.(plot_variables{i})(:,floor(ny/2),iz)},['$$',plot_titles{i},'$$'],'x','');
end
suptitle('last profile at $$y=\pi$$')
print('plot/prof_last','-dpng');
close