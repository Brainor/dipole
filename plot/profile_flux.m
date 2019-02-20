%% flux随时间演化

ix=find(x>.8,1);%最靠近x=1的位置
iz=2;
pei_flux=zeros(size(pei,2)-2,nts);
Gs=zeros(size(pei_flux));
vxs=zeros(size(pei_flux));
Phis=zeros(size(pei_flux));
for nt=1:nts
    load(sprintf('data/dat%4.4d.mat',nt))
    Gs(:,nt)=pei(ix,2:end-1,iz);
    vxs(:,nt)=vex(ix,2:end-1,iz);
    Phis(:,nt)=phi(ix,2:end-1,iz);
end

close all
%% 直接计算flux
pei_flux=Gs.*vxs;
subplot(1,3,1)
draw_plot({1:nts,-mean(pei_flux,1)},'$$\Gamma_G=-Gv_x$$','t');
%% FFT, 1-N
subplot(1,3,2)
for nt=1:nts
    pei_flux(:,nt)=fft(Gs(:,nt))/ny0.*conj(fft(Phis(:,nt))/ny0)*1i.*[0:ny0/2,1-ny0/2:-1]';
end
draw_plot({1:nts,sum(pei_flux,1)},'$$\Gamma_G=\sum_{m=0}^{N-1}im\Phi_mG_m^\star$$','t');
%% FFT, 1-N/2
subplot(1,3,3)
for nt=1:nts
    pei_flux(:,nt)=fft(Gs(:,nt))/ny0.*conj(fft(Phis(:,nt))/ny0)*1i.*(0:ny0-1)';
end
draw_plot({1:nts,sum(real(pei_flux(1:end/2,:)),1)},'$$\Gamma_G=\sum_{m=0}^{N/2}2\Re(im\Phi_mG_m^\star)$$','t');
suptitle('flux of $$G$$ at $$G=0.8$$')