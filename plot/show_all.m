close all
nts_s=[200,200];
%% draw the time evolution
iz=2;%2d
y=-dy:dy:aly;
plot_variables={'phi','pei','deni';...%变量名称
    'wi','G','N';...
    'vex','fluxG','fluxN';...
    'vey','fluxp','fluxn'}';%遍历顺序为先列后行
plot_titles={'\Phi','\widetilde{G}','\widetilde{N}','w','G','N','v_x','\widetilde{G}v_x','\widetilde{N}v_x','v_y','\widetilde{p}_e v_x','\widetilde{n}v_x'};
subplot_grid=num2cell(size(plot_variables'));
close
listings={};
pathname='';
for stage=1:length(nts_s)
    if stage>1
        pathname=[pathname,num2str(stage),'/'];
    end
    filenames=[pathname,'plot/0*.png'];
    if isunix%是UNIX系统
        [~,listing]=unix(['ls -1 ',filenames]);
        listing=split(listing,newline);%为cell
        listing=listing(1:end-1);
    else%是Windows环境
        listing=dir(filenames);
        listing=strcat({listing.folder},'\',{listing.name});
    end
    listings=[listings,listing];
end
makevideo(listings,'time_evolution0');
%% flux随时间演化, movie
plot_variables={'pe0','den0';
    'fluxp','fluxn';
    'p','n'}';%variables name, 遍历顺序为先列后行
plot_titles={'G','N','Gv_x','Nv_x','p_e','\langle n\rangle'};
subplot_grid=num2cell(size(plot_variables'));
v = VideoWriter('flux0','MPEG-4');
v.FrameRate=10;
tic
open(v);
pathname='';
for stage=1:length(nts_s)
    if stage>1
        pathname=[pathname,num2str(stage),'/'];
    end
    for nt=1:nts_s(stage)
        data=load(sprintf([pathname,'data/dat%4.4d.mat'],nt),'pei','vex','deni');
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
        frame = getframe(gcf);
        writeVideo(v,frame);
    end
end
close(v);
time=toc;fprintf('profile一维图: %.2f秒\n',time);