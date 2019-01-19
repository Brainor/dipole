%% p_e与Phi的相位分析
close all
%{
alx=75;nx0=192;
aly=100;ny0=256;
nx=nx0+2;
ny=ny0+2;
%}

ix=find(x>1,1);%最靠近x=1的位置
iz=2;

nk=2;
ky=(nk-1)*2*pi/aly;
st=1;
nts=50;
t=0.001*100*(st:nts);
%{
for nt=1:nts
    if nt<161
        x=80*0.005+x;
    else
        x=100*0.002+x;
    end
    t(nt)=x;
end
%}
petilde=zeros(ny0,1);
vx=petilde;
phitilde=petilde;
pet=zeros(nts-st+1,1);
vxt=pet;
phit=pet;
pek=zeros(ny0,1);
vxk=pek;
phik=pek;
phase=zeros(nts-st+1,1);
sin_phase=phase;
abs_flux=phase;
vy=phase;

for nt=st:nts
    load(sprintf('data/dat%4.4d.mat',nt))
    if (isdeltaf)
        petilde(:)=mean(pei(ix,2:ny0+1,iz));%固定x，取y方向向量
    else
        petilde(:)=pei(ix,2:ny0+1,iz)-mean(pei(ix,2:ny0+1,iz));%固定x，取y方向向量
    end
    vx(:)=vex(ix,2:ny0+1,iz);
    phitilde(:)=phi(ix,2:ny0+1,iz)-mean(phi(ix,2:ny0+1,iz));
    
    pet(nt-st+1)=rms(petilde);%y方向方均根平均
    phit(nt-st+1)=rms(phitilde);
    vxt(nt-st+1)=rms(vx);
    
    pek(:)=fft(petilde(:));
    phik(:)=fft(phitilde(:));
    vxk(:)=fft(vx(:));
    
    phase(nt-st+1)=angle(phik(nk)/pek(nk));
    sin_phase(nt-st+1)=sin(phase(nt-st+1));
    abs_flux(nt-st+1)=ky*abs(pek(nk))*abs(phik(nk))/(ny0*ny0);
    vy(nt-st+1)=mean(mean(vey(:,2:ny0+1,iz),1),2);
    
end
xsep=t(18);%非线性开始
h=subplot(4,1,1);
draw_plot({t(st:nts),phase/pi},'','','\theta_k/\pi','XTickLabel',{},'FontSize',10);
line([xsep,xsep],h.YLim,'Color','k','LineStyle','--');

h=subplot(4,1,2);
draw_plot({t(st:nts),sin_phase},'','','\sin\theta_k','XTickLabel',{},'FontSize',10);
line([xsep,xsep],h.YLim,'Color','k','LineStyle','--');

h=subplot(4,1,3);
draw_plot({t(st:nts),abs_flux},'','','k|p_k||\Phi_k|','XTickLabel',{},'FontSize',10);
line([xsep,xsep],h.YLim,'Color','k','LineStyle','--');

h=subplot(4,1,4);
draw_plot({t(st:nts),vy},'','t(a/c_s)','v_y','FontSize',10);
line([xsep,xsep],h.YLim,'Color','k','LineStyle','--');

suptitle('phase analysis');
print(['plot/phase_flux_phi_nk=',num2str(nk)],'-dpng');
save('data/phase_flux_phi_nk2','t','phase','sin_phase','abs_flux','vy')
close