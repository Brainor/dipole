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
for nt=23:nts%从非线性阶段开始算
    i=i+1;
    load(sprintf('data/dat%4.4d.mat',nt))
    if (isdeltaf)
        pe_mean=pe_mean+mean(pei(:,2:ny0+1,iz)+pe0(:,2:ny0+1,iz),2);
        den_mean=den_mean+mean(deni(:,2:ny0+1,iz)+den0(:,2:ny0+1,iz),2);
        pe_rms=pe_rms+rms(pei(:,2:ny0+1,iz)+pe0(:,2:ny0+1,iz),2);
        den_rms=den_rms+rms(deni(:,2:ny0+1,iz)+den0(:,2:ny0+1,iz),2);
    else
        pe_mean=pe_mean+mean(pei(:,2:ny0+1,iz),2);
        den_mean=den_mean+mean(deni(:,2:ny0+1,iz),2);
        pe_rms=pe_rms+rms(pei(:,2:ny0+1,iz),2);
        den_rms=den_rms+rms(deni(:,2:ny0+1,iz),2);
    end
    phi_mean=phi_mean+mean(phi(:,2:ny0+1,iz),2);
    w_mean=w_mean+mean(wi(:,2:ny0+1,iz),2);
    phi_rms=phi_rms+rms(phi(:,2:ny0+1,iz),2);
    w_rms=w_rms+rms(wi(:,2:ny0+1,iz),2);
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
print('plot/mean_ty','-dpng');
clf;
%% 横坐标为L
for i=1:numel(plot_variables)
    subplot(subplot_grid{:},i);
    draw_plot({1./x,eval([plot_variables{i},'_mean'])},['$$',plot_titles{i},'$$'],'x','');
end
suptitle({'profile with time and $$y$$ average in nonlinear phase','in radial direction'})
print('plot/mean_ty_L','-dpng');
clf;
%% rms平均
for i=1:numel(plot_variables)
    subplot(subplot_grid{:},i);
    draw_plot({x,eval([plot_variables{i},'_rms'])},['$$',plot_titles{i},'$$'],'x','');
end
suptitle('profile with time and $$y$$ rms average in nonlinear phase')
print('plot/mean_ty_rms','-dpng');
close