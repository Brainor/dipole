close all
% global x y
%% draw the time evolution
tic
iz=2;%2d
y=-dy:dy:aly;
plot_variables={'phi','dG','dN';
    'wi','G','N';
    'vey','p','n';
    'flux_G','flux_p','flux_n'};
plot_titles={'\Phi','\widetilde{G}','\widetilde{N}','w','G','N','v_y','p_e','n','Gv_x','p_e v_x','nv_x'};%subplot标题
subplot_grid=num2cell(size(plot_variables'));

parfor nt=1:20
    data=load(sprintf('data/dat%4.4d.mat',nt),'pei','deni','wi','phi','vex','vey');
        data.dG=delt(data.pei);
        data.dN=delt(data.deni);
        data.G=data.pei;
        data.N=data.deni;
    data.n=data.N.*repmat(x',[1,size(deni,2),size(deni,3)]).^4;
    data.p=data.G.*repmat(x',[1,size(deni,2),size(deni,3)]).^(4*gamma);
    data.flux_G=data.G.*data.vex;
    data.flux_p=data.p.*data.vex;
    data.flux_n=data.n.*data.vex;
    h=figure;
    for i=1:numel(plot_variables)
        draw_pcolor([subplot_grid(:)',{i}],{data.(plot_variables{i}),y,x},['$$',plot_titles{i},'$$'],'y','x');
    end
    suptitle(num2str(nt),'FontSize',8);
%     print(sprintf('plot/%4.4d',nt),'-dpng')
savefig(sprintf('plots/%4.4d',nt));
    close(h);
    %     clf;
end
toc
%%
tic
for nt=1:20
%     openfig(sprintf('plots/%4.4d',nt));
uiopen(sprintf('plots/%4.4d.fig',nt),1);
    saveas(gcf,sprintf('plots/%4.4d',nt),'png');
end
toc
%%
close all
tic
iz=2;%2d
y=-dy:dy:aly;
subplot_grid=[4,3];%subplot的前两个参量
parfor nt=1:20
    data=load(sprintf('data/dat%4.4d.mat',nt),'pei','deni','wi','phi','vex','vey');
    data.n=data.deni.*repmat(x',[1,size(deni,2),size(deni,3)]).^4;
    data.p=data.pei.*repmat(x',[1,size(deni,2),size(deni,3)]).^(4*gamma);
    draw_pcolor([4,3,1],{data.phi,y,x},'$$\Phi$$','y','x');
    draw_pcolor([4,3,2],{delt(data.pei),y,x},'$$\widetilde{G}$$','y','x');
    draw_pcolor([4,3,3],{delt(data.deni),y,x},'$$\widetilde{N}$$','y','x');
    draw_pcolor([4,3,4],{data.wi,y,x},'$$w$$','y','x');
    draw_pcolor([4,3,5],{data.pei,y,x},'$$G$$','y','x');
    draw_pcolor([4,3,6],{data.deni,y,x},'$$N$$','y','x');
    draw_pcolor([4,3,7],{data.vey,y,x},'$$v_y$$','y','x');
    draw_pcolor([4,3,8],{data.p,y,x},'$$p_e$$','y','x');
    draw_pcolor([4,3,9],{data.n,y,x},'$$n$$','y','x');
    draw_pcolor([4,3,10],{data.pei.*data.vex,y,x},'$$Gv_x$$','y','x');
    draw_pcolor([4,3,11],{data.p.*data.vex,y,x},'$$p_e v_x$$','y','x');
    draw_pcolor([4,3,12],{data.n.*data.vex,y,x},'$$nv_x$$','y','x');
    suptitle(num2str(nt),'FontSize',8);
    print(sprintf('plot/%4.4d',nt),'-dpng')
    %     clf;
end
toc
%%
% parallel.defaultClusterProfile('local');
% c = parcluster();
job1 = createJob(c);
tic
createTask(job1, @drawtemp, 0, {{1,x,y,gamma},{2,x,y,gamma},{3,x,y,gamma},{4,x,y,gamma},{5,x,y,gamma}});
% % createTask(job1, @drawtemp2, 1, {{1} {2} {3} {4} {5}});
% createTask(job1, @rand, 1, {{3,3} {3,3} {3,3} {3,3} {3,3}});
submit(job1)
wait(job1)
% results=fetchOutputs(job1)
toc
%% 测试plot并行
clc

tic
for i=1:100
    fprintf('%d',i);
    data=load('data/dat0001.mat','pei','deni','wi','vex');
    fields=fieldnames(data);
    for j = 1:numel(fields)
        s=data.(fields{j});
        data.(fields{j})=s(:,12,2);
    end
    h=figure;
    subplot(2,2,1)
    plot(data.deni);
    subplot(2,2,2)
    plot(data.wi);
    subplot(2,2,3)
    plot(data.pei);
    subplot(2,2,4)
    plot(data.vex);
    close(h)
end
fprintf('非并行');
toc
tic
parfor i=1:100
    fprintf('%d',i);
    data=load('data/dat0001.mat','pei','deni','wi','vex');
    fields=fieldnames(data);
    for j = 1:numel(fields)
        s=data.(fields{j});
        data.(fields{j})=s(:,12,2);
    end
    h=figure;
    subplot(2,2,1)
    plot(data.deni);
    subplot(2,2,2)
    plot(data.wi);
    subplot(2,2,3)
    plot(data.pei);
    subplot(2,2,4)
    plot(data.vex);
    close(h)
end
fprintf('并行');
toc
%% 测试pcolor并行
clc
data=load('data/dat0001.mat','pei','deni','wi','vex');
fields=fieldnames(data);
for i = 1:numel(fields)
    s=data.(fields{i});
    data.(fields{i})=s(:,:,2);
end
tic
for i=1:100
    fprintf('%d',i);
    f=figure;
    subplot(2,2,1)
    ph=pcolor(data.deni);shading('interp');ph.ZData=ph.CData;
    subplot(2,2,2)
    ph=pcolor(data.wi);shading('interp');ph.ZData=ph.CData;
    subplot(2,2,3)
    ph=pcolor(data.pei);shading('interp');ph.ZData=ph.CData;
    subplot(2,2,4)
    ph=pcolor(data.vex);shading('interp');ph.ZData=ph.CData;
    close(f);
end
fprintf('\n非并行');
toc
tic
parfor i=1:100
    fprintf('%d',i);
    f=figure;
    subplot(2,2,1)
    ph=pcolor(data.deni);shading('interp');ph.ZData=ph.CData;
    subplot(2,2,2)
    ph=pcolor(data.wi);shading('interp');ph.ZData=ph.CData;
    subplot(2,2,3)
    ph=pcolor(data.pei);shading('interp');ph.ZData=ph.CData;
    subplot(2,2,4)
    ph=pcolor(data.vex);shading('interp');ph.ZData=ph.CData;
    close(f);
end
fprintf('\n并行');
toc

%%
function drawtemp(nt,x,y,gamma)
    fprintf('nt=%d\n',nt);
    data=load(sprintf('data/dat%4.4d.mat',nt),'pei','deni','wi','phi','vex','vey');
    data.n=data.deni.*repmat(x',[1,size(data.deni,2),size(data.deni,3)]).^4;
    data.p=data.pei.*repmat(x',[1,size(data.deni,2),size(data.deni,3)]).^(4*gamma);
    draw_pcolor([4,3,1],{data.phi,y,x},'$$\Phi$$','y','x');
    draw_pcolor([4,3,2],{delt(data.pei),y,x},'$$\widetilde{G}$$','y','x');
    draw_pcolor([4,3,3],{delt(data.deni),y,x},'$$\widetilde{N}$$','y','x');
    draw_pcolor([4,3,4],{data.wi,y,x},'$$w$$','y','x');
    draw_pcolor([4,3,5],{data.pei,y,x},'$$G$$','y','x');
    draw_pcolor([4,3,6],{data.deni,y,x},'$$N$$','y','x');
    draw_pcolor([4,3,7],{data.vey,y,x},'$$v_y$$','y','x');
    draw_pcolor([4,3,8],{data.p,y,x},'$$p_e$$','y','x');
    draw_pcolor([4,3,9],{data.n,y,x},'$$n$$','y','x');
    draw_pcolor([4,3,10],{data.pei.*data.vex,y,x},'$$Gv_x$$','y','x');
    draw_pcolor([4,3,11],{data.p.*data.vex,y,x},'$$p_e v_x$$','y','x');
    draw_pcolor([4,3,12],{data.n.*data.vex,y,x},'$$nv_x$$','y','x');
    suptitle(num2str(nt),'FontSize',8);
    print(sprintf('plot/%4.4d',nt),'-dpng')
end