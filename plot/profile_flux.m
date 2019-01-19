pe_mean=zeros(nx,1);
phi_mean=pe_mean;
den_mean=pe_mean;
vx_mean=pe_mean;
ix=find(x>1,1);%最靠近x=1的位置
iz=2;
for nt=23:nts
    load(sprintf('data/dat%4.4d.mat',nt))
    pe_mean=pe_mean+mean(pei(:,2:ny0+1,iz),2);
    phi_mean=phi_mean+mean(phi(:,2:ny0+1,iz),2);
    den_mean=den_mean+mean(deni(:,2:ny0+1,iz),2);
    vx_mean=vx_mean+mean(vx(:,2:ny0+1,iz),2);
end
pe_mean=pe_mean/nts;
close all
figure
subplot(3,2,1)
draw_plot({x,pe_mean,'b-','Linewidth',1},'$$G$$','','G','XTickLabel',{});
subplot(3,2,3)
draw_plot({x,phi_mean,'b-','Linewidth',1},'$$\Phi$$','','\Phi','XTickLabel',{});
subplot(3,2,5)
draw_plot({x,den_mean,'b-','Linewidth',1},'$$N$$','x','N');
subplot(3,2,2)
draw_plot({x,pe_mean.*x'.^(4*gamma),'b-','Linewidth',1},'$$p_e$$','','p_e','XTickLabel',{});
subplot(3,2,6)
draw_plot({x,den_mean.*x'.^4,'b-','Linewidth',1},'$$\langle n\rangle$$','x','\langle n\rangle');
suptitle('profile with time and $$y$$ average in nonlinear phase','Interpreter','latex')
print('plot/mean_ty','-dpng');