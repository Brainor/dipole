close all
% Phi的频谱图像随时间演化
ix=find(x>1,1);%最靠近x=1的位置
iz=2;
mkdir('temp')
pe0=pe0;
for nt=1:nts
    data=load(sprintf('data/dat%4.4d.mat',nt),'phi','pei');
    phim=amplitude(data.phi(ix,2:end-1,iz));
    if (isdeltaf)
        pem=amplitude(data.pei(ix,2:end-1,iz)+pe0(ix,2:end-1,iz));
    else
        pem=amplitude(data.pei(ix,2:end-1,iz));
    end
    draw_plot({0:ny/2-1,phim,0:ny/2-1,pem},num2str(nt),'m','','YScale','log');
    lgd=legend('$$\Phi$$','$$G$$');lgd.Interpreter='latex';lgd.Box='off';
    print(sprintf('temp/%4.4d',nt),'-dpng');
    clf
end
makevideo('temp/*.png','time_evolution_of_spectrum');
rmdir temp s
%% 非线性时间平均
phim=zeros(1,ny);
pem=phim;
for nt=20:nts
    data=load(sprintf('data/dat%4.4d.mat',nt));
    
    phim=abs(fft(data.phi(ix,:,iz)))/ny+phim;
    if (isdeltaf)
        pem=abs(fft(data.pei(ix,:,iz)+pe0(ix,:,iz)))/ny+pem;
    else
        pem=abs(fft(data.pei(ix,:,iz)))/ny+pem;
    end
end
phim=phim/nts;
pem=pem/nts;
draw_plot({1:floor(ny/2),phim(1:floor(ny/2)),1:floor(ny/2),pem(1:floor(ny/2))},'20-50 slides time average','m','','YScale','log');
lgd=legend('$$\Phi$$','$$G$$');lgd.Interpreter='latex';lgd.Box='off';

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