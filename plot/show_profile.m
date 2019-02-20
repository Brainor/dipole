%% 环向和时间平均profile
plot_variables={'pe','p';...%variables name
    'den','n';...
    'phi','w'}';%遍历顺序为先列后行
plot_titles={'G','p_e','N','\langle n\rangle','\Phi','w'};
subplot_grid=num2cell(size(plot_variables'));
for i=1:numel(plot_variables)%初始化
    eval([plot_variables{i},'_mean=zeros(nx,1);']);
    eval([plot_variables{i},'_rms=zeros(nx,1);']);
end
iz=2;
i=0;
for nt=1:nts%从非线性阶段开始算
    i=i+1;
    load(sprintf('data/dat%4.4d.mat',nt))
    pe_mean=pe_mean+mean(pei(:,2:end-1,iz)+isdeltaf*pe0(:,2:end-1,iz),2);
    den_mean=den_mean+mean(deni(:,2:end-1,iz)+isdeltaf*den0(:,2:end-1,iz),2);
    pe_rms=pe_rms+rms(pei(:,2:end-1,iz)+isdeltaf*pe0(:,2:end-1,iz),2);
    den_rms=den_rms+rms(deni(:,2:end-1,iz)+isdeltaf*den0(:,2:end-1,iz),2);
    phi_mean=phi_mean+mean(phi(:,2:end-1,iz),2);
    w_mean=w_mean+mean(wi(:,2:end-1,iz),2);
    phi_rms=phi_rms+rms(phi(:,2:end-1,iz),2);
    w_rms=w_rms+rms(wi(:,2:end-1,iz),2);
end
pe_mean=pe_mean/i;
p_mean=pe_mean.*x'.^(4*gamma);
den_mean=den_mean/i;
n_mean=den_mean.*x'.^4;
phi_mean=phi_mean/i;
w_mean=w_mean/i;

pe_rms=pe_rms/i;
p_rms=pe_rms.*x'.^(4*gamma);
den_rms=den_rms/i;
n_rms=den_rms.*x'.^4;
phi_rms=phi_rms/i;
w_rms=w_rms/i;
close all
%% 横坐标为x
for i=1:numel(plot_variables)
    subplot(subplot_grid{:},i);
    draw_plot({x,eval([plot_variables{i},'_mean'])},['$$',plot_titles{i},'$$'],'x','');
end
suptitle('profile with time and $$y$$ average in nonlinear phase')
print('plot/profile_x_meanyt','-dpng');
clf;
%% 横坐标为L
for i=1:numel(plot_variables)
    subplot(subplot_grid{:},i);
    draw_plot({1./x,eval([plot_variables{i},'_mean'])},['$$',plot_titles{i},'$$'],'x','');
end
suptitle({'profile with time and $$y$$ average in nonlinear phase','in radial direction'})
print('plot/profile_L_meanyt','-dpng');
clf;
%% rms平均
for i=1:numel(plot_variables)
    subplot(subplot_grid{:},i);
    draw_plot({x,eval([plot_variables{i},'_rms'])},['$$',plot_titles{i},'$$'],'x','');
end
suptitle('profile with time and $$y$$ rms average in nonlinear phase')
print('plot/profile_x_rmsyt','-dpng');
close

%% 环向平均profile
nt_slice=[1,16,30];%初始状态, 线性阶段, 非线性阶段
subplot_grid={3,2};
if (isdeltaf)
    plot_variables={'pei(:,2:end-1,iz)+pe0','p',...%variables name
        'deni(:,2:end-1,iz)+den0','n',...
        'wi','phi'};
else
    plot_variables={'pei','p',...%variables name
        'deni','n',...
        'wi','phi'};
end
plot_titles={'G','p_e','N','\langle n\rangle','w','\Phi'};%variables title
plot_subtitles={'in initial state','in linear state','in nonlinear state'};
figure
for i=1:length(nt_slice)
    clf;
    load(sprintf('data/dat%4.4d.mat',nt_slice(i)))
    if (isdeltaf)
        n=(deni+den0).*repmat(x',[1,size(deni,2),size(deni,3)]).^4;
        p=(pei+pe0).*repmat(x',[1,size(deni,2),size(deni,3)]).^(4*gamma);
    else
        n=deni.*repmat(x',[1,size(deni,2),size(deni,3)]).^4;
        p=pei.*repmat(x',[1,size(deni,2),size(deni,3)]).^(4*gamma);
    end
    for j=1:length(plot_variables)
        subplot(subplot_grid{:},j);
        draw_plot({x,rms(eval([plot_variables{j},'(:,2:end-1,iz)']),2)},['$$',plot_titles{j},'$$'],'x','');
    end
    suptitle({'profile with $$y$$ rms average',plot_subtitles{i}});
    print(['plot/profile_x_rmsy_',num2str(i)],'-dpng');
end
close