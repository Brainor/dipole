% 本征值求解局域不稳定性, 参数由eigen_global给定...
% x, y方向分别为omega_n/omega_d和eta
% 参见Mauel, Bounce-averaged electron fluid equations for interchanged modes
% in a dipole confined plasma
% 和Mauel, Drift orbits in a Symmetric Magnetic Dipole, 2015 12节
%% with N evolution
omegapds=linspace(1,2,500);
% omegapds=linspace(1.6,1.7,200);
etas=linspace(0.1,2,2000);
[X,Y]=meshgrid(etas,omegapds);
rhom=1/25;
gamma=5/3;
m_perp=1;
m=1;

gr=inf(size(X));
for i=1:length(omegapds)
    for j=1:length(etas)
        omegapd=omegapds(i);eta=etas(j);
        A=[0,1,omegapd*1/(1+eta)-1;
            -gamma,2*gamma,omegapd-gamma;
            0,-1/0.6/rhom^2/m_perp^2,0];
        w=eigs(A,1,'largestimag');
        if imag(w)>0
            gr(i,j)=m*max(w);
        end
%         w=imag(eigs(A));
%         w=sort(w);
%         if max(w)>0 
%             gr(i,j)=max(w);
%         end
    end
end
% scan eta, omega_d
figure
h=subplot(1,2,1);
contour(X,Y,imag(gr),'LevelStep',1);
h.XLabel.String='\eta';
h.YLabel.Interpreter='latex';h.YLabel.String='$$\frac{\omega^\star_p}{\omega_d}$$';
h.Title.String='Growth Rate';colorbar;
h=subplot(1,2,2);
contour(X,Y,real(gr),'LevelStep',1);
h.XLabel.String='\eta';
h.YLabel.Interpreter='latex';h.YLabel.String='$$\frac{\omega^\star_p}{\omega_d}$$';
h.Title.String='Frequency';colorbar;
suptitle('extended MHD Mode with heat flux $$\frac{p_e^2}{\tilde{N}}$$','Interpreter','latex');
an=annotation('textbox','LineStyle','none','FontSize',14);
an.Interpreter='latex';
an.String={['$$m=' num2str(m) '$$'],['$$m_\perp=' num2str(m_perp) '$$'],['$$\rho_\star=' num2str(rhom) '$$']};
% scan eta
figure%只画eta
index=find(omegapds>=1.6666,1);
subplot(1,2,1)
h=draw_plot({etas,imag(gr(index,:))},'Growth Rate','\eta','\gamma/m\omega_d');
line([2/3,2/3],h.YLim,'Color','k','LineStyle','--')
subplot(1,2,2)
h=draw_plot({etas,real(gr(index,:))},'Frequency','\eta','\omega/\omega_d');
line([2/3,2/3],h.YLim,'Color','k','LineStyle','--')
line(h.XLim,[0,0],'Color','k','LineStyle','--')
an=annotation('textbox','LineStyle','none','FontSize',14);
an.Interpreter='latex';
an.String={['$$m=' num2str(m) '$$'],['$$m_\perp=' num2str(m_perp) '$$'],['$$\rho_\star=' num2str(rhom) '$$']};
suptitle('extended MHD Mode with heat flux $$\frac{p_e^2}{\tilde{N}}$$','Interpreter','latex');
% close
%% without N evolution
omegapds=linspace(1,2,500);
etas=linspace(0.4,2,200);
[X,Y]=meshgrid(etas,omegapds);
rhom=1/25;
gamma=5/3;
m_perp=1;
m=1;
gr=inf(size(X));
for i=1:length(omegapds)
    for j=1:length(etas)
        omegapd=omegapds(i);eta=etas(j);
        A=[2*gamma,omegapd-gamma;
           -1/0.6/rhom^2/m_perp^2,0];
        w=eigs(A,1,'largestimag');
        if imag(w)>0
            gr(i,j)=m*max(w);
        end
    end
end
figure
h=subplot(1,2,1);
contour(X,Y,imag(gr),'LevelStep',1);
h.XLabel.String='\eta';
h.YLabel.Interpreter='latex';h.YLabel.String='$$\frac{\omega^\star_p}{\omega_d}$$';
h.Title.String='Growth Rate';colorbar;
h=subplot(1,2,2);
contour(X,Y,real(gr),'LineWidth',2);
h.XLabel.String='\eta';
h.YLabel.Interpreter='latex';h.YLabel.String='$$\frac{\omega^\star_p}{\omega_d}$$';
h.Title.String='Frequency';colorbar;
suptitle('ideal MHD Mode with heat flux $$p_e^2$$','Interpreter','latex');
close
%% with N, scan m,omegapds
omegapds=linspace(1,2,500);
omegapd0=[2,1.9,1.8,1.7,1.6,1.5];%取特定omegapd
iy=zeros(size(omegapd0));
eta=0;
ms=1:80;
[X,~]=meshgrid(ms,omegapds);
rhom=1/25;
gamma=5/3;
m=1;

gr=inf(size(X));
for i=1:length(omegapds)
    for j=1:length(ms)
        omegapd=omegapds(i);m=ms(j);
        A=[0,1,omegapd*1/(1+eta)-1;
            -gamma,2*gamma,omegapd-gamma;
            0,-1/0.6/rhom^2/m^2,0];
        w=eigs(A,1,'largestimag');
        if imag(w)>0
            gr(i,j)=m*max(w);
        end
    end
end
for i=1:length(iy)
    iy(i)=find(omegapds>=omegapd0(i),1);
end
subplot(2,2,1)
draw_plot({ms*rhom,imag(gr(iy,:)),'.','Linewidth',1},...
    'growth rate','k_y \rho_i','\gamma/\omega_d^{-1}');
lgd=legend(cellfun(@num2str,num2cell(omegapd0),'UniformOutput',false),'Location','best');
lgd.Title.String='$$\frac{\omega_p^\star}{\omega_d}$$';lgd.Title.Interpreter='latex';
subplot(2,2,2)
draw_plot({ms*rhom,real(gr(iy,:)),'.','Linewidth',1},...
    'frequency','k_y \rho_i','\omega/\omega_d^{-1}');
lgd=legend(cellfun(@num2str,num2cell(omegapd0),'UniformOutput',false),'Location','best');
lgd.Title.String='$$\frac{\omega_p^\star}{\omega_d}$$';lgd.Title.Interpreter='latex';
subplot(2,2,3)
draw_plot({ms,imag(gr(iy,:)),'.','Linewidth',1},...
    'growth rate','m','\gamma/\omega_d^{-1}');
lgd=legend(cellfun(@num2str,num2cell(omegapd0),'UniformOutput',false),'Location','best');
lgd.Title.String='$$\frac{\omega_p^\star}{\omega_d}$$';lgd.Title.Interpreter='latex';
subplot(2,2,4)
draw_plot({ms,real(gr(iy,:)),'.','Linewidth',1},...
    'frequency','m','\omega/\omega_d^{-1}');
lgd=legend(cellfun(@num2str,num2cell(omegapd0),'UniformOutput',false),'Location','best');
lgd.Title.String='$$\frac{\omega_p^\star}{\omega_d}$$';lgd.Title.Interpreter='latex';
suptitle({'dispersion relationship';['$$\eta=',num2str(eta),'$$, $$\rho_\star=' num2str(rhom) '$$']});
print(['plot/growthrate_eta=',num2str(eta),'.png'],'-dpng');
%% with N, 计算flux
% Gamma_N=real(i*m*Phi*N')
% Gamma_p=real(i*m*Phi*p*)
omegapds=linspace(1,2,500);omegapds=linspace(1.2,3,1000);
etas=linspace(0.4,2,200);etas=linspace(0,6,400);
[X,Y]=meshgrid(etas,omegapds);
rhom=1/25;
gamma=5/3;
m_perp=1;
m=1;

gr=nan(size(X));
vectors=nan(size(X,1),size(X,2),3);%N, p, Phi
for i=1:length(omegapds)
    for j=1:length(etas)
        omegapd=omegapds(i);eta=etas(j);
        A=[0,1,omegapd*1/(1+eta)-1;
            -gamma,2*gamma,omegapd-gamma;
            0,-1/0.6/rhom^2/m_perp^2,0];
        [V,w]=eigs(A,1,'largestimag');%分别为特征向量和特征值
        if imag(w)>0
            gr(i,j)=m*max(w);
            vectors(i,j,:)=V;
        end
%         w=imag(eigs(A));
%         w=sort(w);
%         if max(w)>0 
%             gr(i,j)=max(w);
%         end
    end
end
h=subplot(2,2,1);
ph=pcolor(X,Y,real(1i*m*vectors(:,:,3).*conj(vectors(:,:,1))));shading interp;ph.ZData=ph.CData;
h.XLabel.String='\eta';
h.YLabel.Interpreter='latex';h.YLabel.String='$$\frac{\omega^\star_p}{\omega_d}$$';
h.Title.Interpreter='latex';h.Title.String={'Particle flux';'$\Gamma_N=\Re\{im\Phi_mP_m^\star\}$'};colorbar;
h=subplot(2,2,2);
ph=pcolor(X,Y,real(1i*m*vectors(:,:,3).*conj(vectors(:,:,2))));shading interp;ph.ZData=ph.CData;
h.XLabel.String='\eta';
h.YLabel.Interpreter='latex';h.YLabel.String='$$\frac{\omega^\star_p}{\omega_d}$$';
h.Title.Interpreter='latex';h.Title.String={'Radial pressure flux';'$\Gamma_p=\Re\{im\Phi_mp_m^\star\}$'};colorbar;
h=subplot(2,2,3);
test=angle(vectors(:,:,3)./vectors(:,:,1));
test(test>0)=test(test>0)-2*pi;
ph=pcolor(X,Y,test);shading interp;ph.ZData=ph.CData;
h.XLabel.String='\eta';
h.YLabel.Interpreter='latex';h.YLabel.String='$$\frac{\omega^\star_p}{\omega_d}$$';
h.Title.String='Phase lag between $\Phi_m$, $N_m$';h.Title.Interpreter='latex';colorbar;
h=subplot(2,2,4);
ph=pcolor(X,Y,angle(vectors(:,:,3)./vectors(:,:,2)));shading interp;ph.ZData=ph.CData;
h.XLabel.String='\eta';
h.YLabel.Interpreter='latex';h.YLabel.String='$$\frac{\omega^\star_p}{\omega_d}$$';
h.Title.String='Phase lag between $\Phi_m$, $p_m$';h.Title.Interpreter='latex';colorbar;
% annotation('textbox','LineStyle','none','Interpreter','latex','FontSize',12,'String',{['$$m=' num2str(m) '$$'],['$$m_\perp=' num2str(m_perp) '$$'],['$$\rho_\star=' num2str(rhom) '$$']});
