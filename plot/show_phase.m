%% p_e与Phi的相位分析
close all

ix=find(x>1,1);%最靠近x=1的位置
iz=2;
m=2;%模数
st=36;%非线性的位置


plot_variables={'phase_G','phase_N';
    'phase_G_sin','phase_N_sin';
    'flux_G_abs','flux_N_abs';
    'flux_G','flux_N';
    'vy','vy'}';
plot_titles={'\theta_\Phi-\theta_G','\theta_\Phi-\theta_N',...
    '\sin(\theta_\Phi-\theta_G)','\sin(\theta_\Phi-\theta_N)',...
    'm|G_m||\Phi_m|','m|N_m||\Phi_m|',...
    '2\Re(im\Phi_m\delta G_m^*)','2\Re(im\Phi_m\delta N_m^*)',...
    'v_y','v_y'};

subplot_grid=num2cell(size(plot_variables'));

t=ntp*tau*(1:nts);
baseline=t(st);

data2.phase_G=zeros(nts,1);
data2.phase_N=data2.phase_G;
data2.flux_G_abs=data2.phase_G;%|Gamma|=|mPhi*N|
data2.flux_G=data2.phase_G;%Gamma=2Re(imPhi*N')
data2.flux_N_abs=data2.phase_G;data2.flux_N=data2.phase_G;
data2.phase_G_sin=data2.phase_G;data2.phase_N_sin=data2.phase_G;
data2.vy=data2.phase_G;


for nt=1:nts
    data=load(sprintf('data/dat%4.4d.mat',nt),'pei','phi','phi','vey','deni');
    if (isdeltaf)
        data.dG=data.pei(ix,2:ny0+1,iz);%delta G
        data.dN=data.deni(ix,2:ny0+1,iz);%delta N
    else
        data.dG=delt(data.pei(ix,2:ny0+1,iz),2);
        data.dN=delt(data.deni(ix,2:ny0+1,iz),2);
    end
    data.dPhi=delt(data.phi(ix,2:ny0+1,iz),2);
    
    data.Gm=fft(data.dG)/ny0;
    data.Nm=fft(data.dN)/ny0;
    data.Phim=fft(data.dPhi)/ny0;
    
    data2.phase_G(nt)=angle(data.Phim(m+1)/data.Gm(m+1));
    data2.phase_N(nt)=angle(data.Phim(m+1)/data.Nm(m+1));
    data2.vy(nt)=mean(data.vey(:,2:ny0+1,iz),'all');
    data2.flux_G_abs(nt)=m*abs(data.Gm(m+1)*data.Phim(m+1));
    data2.flux_G(nt)=2*real(1i*m*data.Phim(m+1).*conj(data.Gm(m+1)));
    data2.flux_N_abs(nt)=m*abs(data.Nm(m+1)*data.Phim(m+1));
    data2.flux_N(nt)=2*real(1i*m*data.Phim(m+1).*conj(data.Nm(m+1))); 
end
data2.phase_G_sin=sin(data2.phase_G);data2.phase_N_sin=sin(data2.phase_N);
%%
for i=1:numel(plot_variables)
    subplot(subplot_grid{:},i);
    h=draw_plot({t,data2.(plot_variables{i})},['$$',plot_titles{i},'$$']);
    line(baseline*ones(1,2),h.YLim,'LineStyle','--','Color','k');
    if mod(i,5)~=0
        h.XTickLabel={};
    end
end

suptitle({'phase analysis';['$$m=',num2str(m),'$$']});
print(['plot/phase_m=',num2str(m)],'-dpng');