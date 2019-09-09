close all
% Phi的频谱图像随时间演化
ix1=find(x>1,1);%最靠近x=1的位置
ix2=find(x>0.7,1);
iz=2;
mkdir('temp')
pe0=pe0;
den0=den0;
for nt=150:nts
    data=load(sprintf('data/dat%4.4d.mat',nt),'phi','pei','deni');
    phim1=amplitude(data.phi(ix1,2:end-1,iz));
    phim2=amplitude(data.phi(ix2,2:end-1,iz));
    if (isdeltaf)
        denm1=amplitude(data.deni(ix,2:end-1,iz)+den0(ix,2:end-1,iz));
        pem1=amplitude(data.pei(ix,2:end-1,iz)+pe0(ix,2:end-1,iz));
    else
        denm1=amplitude(data.deni(ix1,2:end-1,iz));
        pem1=amplitude(data.pei(ix1,2:end-1,iz));
        denm2=amplitude(data.deni(ix2,2:end-1,iz));
        pem2=amplitude(data.pei(ix2,2:end-1,iz));
    end
    %     draw_plot({0:ny/2-1,[phim1;pem1;denm1;phim2;pem2;denm2]},num2str(nt),'m','','YScale','log');
    draw_plot({0:40,[phim1(1:41);pem1(1:41);denm1(1:41);phim2(1:41);pem2(1:41);denm2(1:41)]},num2str(nt),'m','','YScale','log');
    lgd=legend('$$\Phi(x=1)$$','$$G(x=1)$$','$$N(x=1)$$','$$\Phi(x=0.5)$$','$$G(x=0.5)$$','$$N(x=0.5)$$');lgd.Interpreter='latex';lgd.Box='off';
    print(sprintf('temp/%4.4d',nt),'-dpng');
    clf
end
makevideo('temp/*.png','time_evolution_of_spectrums');
rmdir temp s
%% 非线性时间平均
nt=50;
data=load(sprintf('data/dat%4.4d.mat',nt),'phi','pei','deni');
phim1=amplitude(data.phi(ix1,2:end-1,iz));
phim2=amplitude(data.phi(ix2,2:end-1,iz));
denm1=amplitude(data.deni(ix1,2:end-1,iz));
pem1=amplitude(data.pei(ix1,2:end-1,iz));
denm2=amplitude(data.deni(ix2,2:end-1,iz));
pem2=amplitude(data.pei(ix2,2:end-1,iz));

denm1(1)=inf;
denm2(1)=inf;
subplot(1,2,1)
draw_plot({0:59,[phim1(1:60);denm1(1:60);phim2(1:60);denm2(1:60)]},'50','m','','YScale','log','Box','on');

phim1=zeros(1,floor(ny/2));
pem1=phim1;denm1=phim1;
phim2=phim1;
pem2=phim1;denm2=phim1;
for nt=150:nts
    data=load(sprintf('data/dat%4.4d.mat',nt));
    
    phim1=amplitude(data.phi(ix1,2:end-1,iz))+phim1;
    phim2=amplitude(data.phi(ix2,2:end-1,iz))+phim2;
    if (isdeltaf)
        pem=abs(fft(data.pei(ix,:,iz)+pe0(ix,:,iz)))/ny+pem;
        denm1=abs(fft(data.deni(ix,:,iz)+den0(ix,:,iz)))/ny+denm1;
    else
        pem1=amplitude(data.pei(ix1,2:end-1,iz))+pem1;
        denm1=amplitude(data.deni(ix1,2:end-1,iz))+denm1;
        pem2=amplitude(data.pei(ix2,2:end-1,iz))+pem2;
        denm2=amplitude(data.deni(ix2,2:end-1,iz))+denm2;
    end
end
phim1=phim1/(nts-150+1);
pem1=pem1/(nts-150+1);
denm1=denm1/(nts-150+1);
phim2=phim2/(nts-150+1);
pem2=pem2/(nts-150+1);
denm2=denm2/(nts-150+1);
denm1(1)=inf;
denm2(1)=inf;
subplot(1,2,2)
draw_plot({1:floor(ny/2),[phim1(1:floor(ny/2));pem1(1:floor(ny/2));denm1(1:floor(ny/2));phim2(1:floor(ny/2));pem2(1:floor(ny/2));denm2(1:floor(ny/2))]},['100-',num2str(nts),' slides time average'],'m','','YScale','log');
lgd=legend('$$\Phi(x=1)$$','$$G(x=1)$$','$$N(x=1)$$','$$\Phi(x=0.5)$$','$$G(x=0.5)$$','$$N(x=0.5)$$');lgd.Interpreter='latex';lgd.Box='off';
draw_plot({0:59,[phim1(1:60);denm1(1:60);phim2(1:60);denm2(1:60)]},['150-',num2str(nts),' slides time average'],'m','','YScale','log','Box','on');
lgd=legend('$$\Phi(x=1)$$','$$N(x=1)$$','$$\Phi(x=0.7)$$','$$N(x=0.7)$$');lgd.Interpreter='latex';lgd.Box='off';

print(['plot/spec_x=',sprintf('%.2f',x(ix)),'.png'],'-dpng');
function vm=amplitude(v)
    % AMPLITUDE calculate the amplitude of phase
    % v is the real vector
    % vm has length of floor(n/2)+1, with phase number 0:floor(n/2)+1
    n=length(v);
    vm=abs(fft(v))/n;
    vm(2:floor((n+1)/2))=2*vm(2:floor((n+1)/2));
    vm(floor(n/2)+2:end)=[];
end