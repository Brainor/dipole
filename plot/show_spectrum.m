close all
% Phi的频谱图像随时间演化
ix=find(x>1,1);%最靠近x=1的位置
% ix=find(x>0.5,1);
iz=2;
mkdir('temp')
pe0=pe0;
den0=den0;
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
phim=zeros(1,ny);
pem=phim;denm=phim;
for nt=1:nts
    data=load(sprintf('data/dat%4.4d.mat',nt));
    
    phim=abs(fft(data.phi(ix,:,iz)))/ny+phim;
    if (isdeltaf)
        pem=abs(fft(data.pei(ix,:,iz)+pe0(ix,:,iz)))/ny+pem;
        denm=abs(fft(data.deni(ix,:,iz)+den0(ix,:,iz)))/ny+denm;
    else
        pem=abs(fft(data.pei(ix,:,iz)))/ny+pem;
        denm=abs(fft(data.deni(ix,:,iz)))/ny+denm;
    end
end
phim=phim/nts;
pem=pem/nts;
denm=denm/nts;
draw_plot({1:floor(ny/2),[phim(1:floor(ny/2));pem(1:floor(ny/2));denm(1:floor(ny/2))]},['1-',num2str(nts),' slides time average'],'m','','YScale','log');
lgd=legend('$$\Phi$$','$$G$$','$$N$$');lgd.Interpreter='latex';lgd.Box='off';

print(['plot/spec_x=',sprintf('%.2f',x(ix)),'.png'],'-dpng');
close
function vm=amplitude(v)
% AMPLITUDE calculate the amplitude of phase
% v is the real vector
% vm has length of floor(n/2)+1, with phase number 0:floor(n/2)+1
n=length(v);
vm=abs(fft(v))/n;
vm(2:floor((n+1)/2))=2*vm(2:floor((n+1)/2));
vm(floor(n/2)+2:end)=[];
end