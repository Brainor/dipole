close all
% Phi的频谱图像随时间演化
ix=find(x>1,1);%最靠近x=1的位置
% ix=find(x>0.5,1);
iz=2;
mkdir('temp')
pe0=pe0;
den0=den0;
close
for nt=1:nts
    data=load(sprintf('data/dat%4.4d.mat',nt),'phi','pei','deni');
    phim=amplitude(data.phi(ix,2:end-1,iz));
    if (isdeltaf)
        denm=amplitude(data.deni(ix,2:end-1,iz)+den0(ix,2:end-1,iz));
        pem=amplitude(data.pei(ix,2:end-1,iz)+pe0(ix,2:end-1,iz));
    else
        denm=amplitude(data.deni(ix,2:end-1,iz));
        pem=amplitude(data.pei(ix,2:end-1,iz));
    end
    draw_plot({0:ny/2-1,[phim;pem;denm]},num2str(nt),'m','','YScale','log');
    lgd=legend('$$\Phi$$','$$G$$','$$N$$');lgd.Interpreter='latex';lgd.Box='off';
    print(sprintf('temp/%4.4d',nt),'-dpng');
    clf
end
makevideo('temp/*.png','spectrum_t');
rmdir temp s
%% 非线性时间平均
phim=zeros(1,ny0);
pem=phim;denm=phim;
for nt=1:nts
    data=load(sprintf('data/dat%4.4d.mat',nt));
    
    phim=abs(fft(data.phi(ix,2:end-1,iz)))/ny+phim;
    if (isdeltaf)
        pem=abs(fft(data.pei(ix,2:end-1,iz)+pe0(ix,:,iz)))/ny+pem;
        denm=abs(fft(data.deni(ix,2:end-1,iz)+den0(ix,:,iz)))/ny+denm;
    else
        pem=abs(fft(data.pei(ix,2:end-1,iz)))/ny+pem;
        denm=abs(fft(data.deni(ix,2:end-1,iz)))/ny+denm;
    end
end
phim=phim/nts;
pem=pem/nts;
denm=denm/nts;
draw_plot({1:floor(ny/2),[phim(1:floor(ny/2));pem(1:floor(ny/2));denm(1:floor(ny/2))]},['1-',num2str(nts),' slides time average'],'m','','YScale','log');
lgd=legend('$$\Phi$$','$$G$$','$$N$$');lgd.Interpreter='latex';lgd.Box='off';

print(['plot/spec_x=',sprintf('%.2f',x(ix)),'.png'],'-dpng');
%% 模结构
ms=[1,2,4];%取这几个模
nt=40;%选取时间点
data=load(sprintf('data/dat%4.4d.mat',nt));
phim=data.phi(:,2:end-1,2);
phim=real(fft(phim,[],2))/ny;
phim(2:end)=2*phim(2:end);
phi_m=phim(:,ms+1);%要做图的phim
% subplot(1,2,1)
% draw_plot({x,phi_m});
% legend(cellfun(@num2str,num2cell(ms),'UniformOutput',false),'Location','best','Box','off');
% subplot(1,2,2)%与eigen作比较
phim_eigen([1,3],:)=-1*phim_eigen([1,3],:);

for i=1:length(ms)%归一化
    phim_eigen(i,:)=max(abs(phi_m(:,i)))/max(abs(phim_eigen(i,:)))*phim_eigen(i,:);
end
close
hold on
h=draw_plot({x(2:end-1),phi_m(2:end-1,:)},'','','','ColorOrderIndex',1);
legend(h,cellfun(@num2str,num2cell(ms),'UniformOutput',false),'Location','northeast','Box','off');
h2=copyobj(h,gcf);%delete(get(h2,'Children'));
hold on;draw_plot({xx(2:end-1),real(phim_eigen),'--','Parent',h2},['nt=',num2str(nt)],'x');
axis(h2,'off');
legend([h2.Children(length(ms)+1),h2.Children(1)],{'ivp','eigen'},'Location','northwest','Box','off');
hold off
print('plot/eigen','-dpng');
%% 时空频谱
x_s=[0.7,1];
ix=arrayfun(@(x_s)find(x>x_s,1),x_s);
data=struct;data2=struct;data3=struct;data4=struct;
for nt=150:nts
    data=load(sprintf('data/dat%4.4d.mat',nt),'phi','pei','deni');
    data2.phi(1:length(x_s),1:ny0,nt-150+1)=data.phi(ix,2:end-1,iz);
    data2.pei(1:length(x_s),1:ny0,nt-150+1)=data.pei(ix,2:end-1,iz);
    data2.deni(1:length(x_s),1:ny0,nt-150+1)=data.deni(ix,2:end-1,iz);
end
% 处理空域
for i=1:length(x_s)
    for j=1:nts-150+1
        data3.phi(i,:,j)=amplitude(data2.phi(i,:,j));
        data3.pei(i,:,j)=amplitude(data2.pei(i,:,j));
        data3.deni(i,:,j)=amplitude(data2.deni(i,:,j));
    end
end
%处理时域
for i=1:length(x_s)
    for j=1:ny0/2+1
        data4.phi(i,j,:)=abs(ifft(data3.phi(i,j,:)));
        data4.pei(i,j,:)=abs(ifft(data3.pei(i,j,:)));
        data4.deni(i,j,:)=abs(ifft(data3.deni(i,j,:)));
    end
end
h=subplot(1,2,1);
draw_pcolor(h,{squeeze(data4.phi(1,:,:)),1:51,1:ny0/2+1},'','m','f','FontSize',12);
h=subplot(1,2,2);
h=draw_pcolor(h,{squeeze(data4.phi(2,:,:)),1:51,1:ny0/2+1},'','m','f','FontSize',12);

function vm=amplitude(v)
% AMPLITUDE calculate the amplitude of phase
% v is the real vector
% vm has length of floor(n/2)+1, with phase number 0:floor(n/2)+1
n=length(v);
vm=abs(fft(v))/n;
vm(2:floor((n+1)/2))=2*vm(2:floor((n+1)/2));
vm(floor(n/2)+2:end)=[];
end