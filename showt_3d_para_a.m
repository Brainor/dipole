close all

dx=alx/(nx-1);
x=dx*(0:nx-1);

dy=aly/(ny-1);
y=dy*(0:ny-1);

[X,Y]=meshgrid(y,x);

lx=[nx2/nx0*alx,nx2/nx0*alx];
ly=[0,aly];

x_sol=nx2/nx0*alx;

datlist=dir('dat*.mat');
nts=length(datlist);

izOut=6;
izIn=14;

iy=10;
ix=10;
% parpool('local',4);

parfor nt=1:nts
    data=load(['dat',sprintf('%4.4d',nt)]);

    % plot x-y outside
    figure;
    subplot(421)
    pcolor(X,Y,data.phi(:,:,izOut)); colorbar; shading interp; colormap jet;
    title('$\phi$','Interpreter','latex');
    ylabel('$x/\rho_s$','Interpreter','latex');
    xlabel('$y/\rho_s$','Interpreter','latex');
    hold on;
    plot(ly,lx,'k-','LineWidth',1);

    subplot(422)
    pcolor(X,Y,-data.wi(:,:,izOut)); colorbar; shading interp; colormap jet;
    title('$w$','Interpreter','latex');
    ylabel('$x/\rho_s$','Interpreter','latex');
    xlabel('$y/\rho_s$','Interpreter','latex');
    hold on;
    plot(ly,lx,'k-','LineWidth',1);

    subplot(423)
    pcolor(X,Y,data.pei(:,:,izOut)); colorbar;shading interp; colormap jet;
    title('$p_e$','Interpreter','latex');
    ylabel('$x/\rho_s$','Interpreter','latex');
    xlabel('$y/\rho_s$','Interpreter','latex');
    hold on;
    plot(ly,lx,'k-','LineWidth',1);

    subplot(424)
    a=sp0(data.pei);
    pcolor(X,Y,a(:,:,izOut));colorbar; shading interp; colormap jet;
    title('${\tilde p}_e$','Interpreter','latex');
    xlabel('$y/\rho_s$','Interpreter','latex');
    ylabel('$x/\rho_s$','Interpreter','latex');
    hold on;
    plot(ly,lx,'k-','LineWidth',1);

    subplot(425)
    pcolor(X,Y,data.vex(:,:,izOut)); colorbar; shading interp; colormap jet;
    title('$v_x$','Interpreter','latex');
    ylabel('$x/\rho_s$','Interpreter','latex');
    xlabel('$y/\rho_s$','Interpreter','latex');
    hold on;
    plot(ly,lx,'k-','LineWidth',1);

    subplot(426)
    pcolor(X,Y,data.vey(:,:,izOut)); colorbar; shading interp; colormap jet;
    title('$v_y$','Interpreter','latex');
    xlabel('$y/\rho_s$','Interpreter','latex');
    ylabel('$x/\rho_s$','Interpreter','latex');
    hold on;
    plot(ly,lx,'k-','LineWidth',1);

    subplot(427)
    pcolor(X,Y,data.vii(:,:,izOut)); colorbar;shading interp; colormap jet;
    title('$u_i$','Interpreter','latex');
    xlabel('$y/\rho_s$','Interpreter','latex');
    ylabel('$x/\rho_s$','Interpreter','latex');

    subplot(428)
    pcolor(X,Y,data.jz(:,:,izOut)); colorbar;shading interp; colormap jet;
    title('$j_z$','Interpreter','latex');
    xlabel('$y/\rho_s$','Interpreter','latex');
    ylabel('$x/\rho_s$','Interpreter','latex');

    filename=sprintf('%2.2dz%4.4d',izOut,nt);
    % saveas(gcf,filename);     % for Matlab r2014a save as .fig
    saveas(gcf,[filename,'.png']);   % for Matlab r2014b save as .png
    close;


    % plot x-y inside
    figure;
    subplot(421)
    pcolor(X,Y,data.phi(:,:,izIn)); colorbar; shading interp; colormap jet;
    title('$\phi$','Interpreter','latex');
    ylabel('$x/\rho_s$','Interpreter','latex');
    xlabel('$y/\rho_s$','Interpreter','latex');
    hold on;
    plot(ly,lx,'k-','LineWidth',1);

    subplot(422)
    pcolor(X,Y,-data.wi(:,:,izIn)); colorbar; shading interp; colormap jet;
    title('$w$','Interpreter','latex');
    ylabel('$x/\rho_s$','Interpreter','latex');
    xlabel('$y/\rho_s$','Interpreter','latex');
    hold on;
    plot(ly,lx,'k-','LineWidth',1);

    subplot(423)
    pcolor(X,Y,data.pei(:,:,izIn)); colorbar;shading interp; colormap jet;
    title('$p_e$','Interpreter','latex');
    ylabel('$x/\rho_s$','Interpreter','latex');
    xlabel('$y/\rho_s$','Interpreter','latex');
    hold on;
    plot(ly,lx,'k-','LineWidth',1);

    subplot(424)
    a=sp0(data.pei);
    pcolor(X,Y,a(:,:,izIn));colorbar; shading interp; colormap jet;
    title('${\tilde p}_e$','Interpreter','latex');
    xlabel('$y/\rho_s$','Interpreter','latex');
    ylabel('$x/\rho_s$','Interpreter','latex');
    hold on;
    plot(ly,lx,'k-','LineWidth',1);

    subplot(425)
    pcolor(X,Y,data.vex(:,:,izIn)); colorbar; shading interp; colormap jet;
    title('$v_x$','Interpreter','latex');
    ylabel('$x/\rho_s$','Interpreter','latex');
    xlabel('$y/\rho_s$','Interpreter','latex');
    hold on;
    plot(ly,lx,'k-','LineWidth',1);

    subplot(426)
    pcolor(X,Y,data.vey(:,:,izIn)); colorbar; shading interp; colormap jet;
    title('$v_y$','Interpreter','latex');
    xlabel('$y/\rho_s$','Interpreter','latex');
    ylabel('$x/\rho_s$','Interpreter','latex');
    hold on;
    plot(ly,lx,'k-','LineWidth',1);

    subplot(427)
    pcolor(X,Y,data.vii(:,:,izIn)); colorbar;shading interp; colormap jet;
    title('$u_i$','Interpreter','latex');
    xlabel('$y/\rho_s$','Interpreter','latex');
    ylabel('$x/\rho_s$','Interpreter','latex');
    hold on;
    plot(ly,lx,'k-','LineWidth',1);

    subplot(428)
    pcolor(X,Y,data.jz(:,:,izIn)); colorbar;shading interp;colormap jet;
    title('$j_z$','Interpreter','latex');
    ylabel('$y/\rho_s$','Interpreter','latex');
    xlabel('$x/\rho_s$','Interpreter','latex');
    hold on;
    plot(ly,lx,'k-','LineWidth',1);

    filename=sprintf('%2.2dz%4.4d',izIn,nt);
    % saveas(gcf,filename);     % for Matlab r2014a save as .fig
    saveas(gcf,[filename,'.png']);   % for Matlab r2014b save as .png
    close;


    % plot x-z
    figure;
    subplot(421)
    pcolor(squeeze(data.phi(:,iy,:))); colorbar; shading interp; colormap jet;
    title('$\phi$','Interpreter','latex');
    ylabel('x');
    xlabel('z');
    hold on;

    subplot(422)
    pcolor(squeeze(-data.wi(:,iy,:))); colorbar; shading interp; colormap jet;
    title('$w$','Interpreter','latex');
    ylabel('x');
    xlabel('z');
    hold on;

    subplot(423)
    pcolor(squeeze(data.pei(:,iy,:))); colorbar;shading interp; colormap jet;
    title('$p_e$','Interpreter','latex');
    ylabel('x');
    xlabel('z');
    hold on;

    subplot(424)
    a=sp0(data.pei);
    pcolor(squeeze(a(:,iy,:)));colorbar; shading interp; colormap jet;
    title('${\tilde p}_e$','Interpreter','latex');
    ylabel('x');
    xlabel('z');
    hold on;

    subplot(425)
    pcolor(squeeze(data.vex(:,iy,:))); colorbar; shading interp; colormap jet;
    title('$v_x$','Interpreter','latex');
    ylabel('x');
    xlabel('z');
    hold on;

    subplot(426)
    pcolor(squeeze(data.vey(:,iy,:))); colorbar; shading interp; colormap jet;
    title('$v_y$','Interpreter','latex');
    ylabel('x');
    xlabel('z');
    hold on;

    subplot(427)
    pcolor(squeeze(data.vii(:,iy,:))); colorbar;shading interp; colormap jet;
    title('$u_i$','Interpreter','latex');
    ylabel('x');
    xlabel('z');
    hold on;

    subplot(428)
    pcolor(squeeze(data.jz(:,iy,:))); colorbar;shading interp; colormap jet;
    title('$j_z$','Interpreter','latex');
    ylabel('x');
    xlabel('z');

    filename=sprintf('%2.2dy%4.4d',iy,nt);
    % saveas(gcf,filename);     % for Matlab r2014a save as .fig
    saveas(gcf,[filename,'.png']);   % for Matlab r2014b save as .png
    close;


   % plot y-z
    figure;
    subplot(421)
    pcolor(squeeze(data.phi(ix,:,:))); colorbar; shading interp; colormap jet;
    title('$\phi$','Interpreter','latex');
    ylabel('y');
    xlabel('z');
    hold on;

    subplot(422)
    pcolor(squeeze(-data.wi(ix,:,:))); colorbar; shading interp; colormap jet;
    title('$w$','Interpreter','latex');
    ylabel('y');
    xlabel('z');
    hold on;

    subplot(423)
    pcolor(squeeze(data.pei(ix,:,:))); colorbar;shading interp; colormap jet;
    title('$p_e$','Interpreter','latex');
    ylabel('y');
    xlabel('z');
    hold on;

    subplot(424)
    a=sp0(data.pei);
    pcolor(squeeze(a(ix,:,:)));colorbar; shading interp; colormap jet;
    title('${\tilde p}_e$','Interpreter','latex');
    ylabel('y');
    xlabel('z');
    hold on;

    subplot(425)
    pcolor(squeeze(data.vex(ix,:,:))); colorbar; shading interp; colormap jet;
    title('$v_x$','Interpreter','latex');
    ylabel('y');
    xlabel('z');
    hold on;

    subplot(426)
    pcolor(squeeze(data.vey(ix,:,:))); colorbar; shading interp; colormap jet;
    title('$v_y$','Interpreter','latex');
    ylabel('y');
    xlabel('z');
    hold on;

    subplot(427)
    pcolor(squeeze(data.vii(ix,:,:))); colorbar;shading interp; colormap jet;
    title('$u_i$','Interpreter','latex');
    ylabel('y');
    xlabel('z');
    hold on;

    subplot(428)
    pcolor(squeeze(data.jz(:,iy,:))); colorbar;shading interp;colormap jet;
    title('$j_z$','Interpreter','latex');
    ylabel('y');
    xlabel('z');
    hold on;

    filename=sprintf('%2.2dx%4.4d',ix,nt);
    % saveas(gcf,filename);     % for Matlab r2014a save as .fig
    saveas(gcf,[filename,'.png']);   % for Matlab r2014b save as .png
    close;


    % % plot profile
    % figure;
    % set(gca,'FontSize',14);
    % subplot(311)
    % plot(x,mean(data.phi(:,:,izOut),2),'r--','Linewidth',1,'Marker','o');
    % hold on;
    % plot(x,mean(data.phi(:,:,izIn),2),'b--','Linewidth',1,'Marker','o');
    % title('Protential Profile');
    % xlabel('x');
    % ylabel('phi');
    % xlim([0,alx]);
    % grid on;
    % set(gca,'xTick',x_sol);
    % set(gca,'LineWidth',1);
    % set(gca,'yGrid','off');

    % subplot(312)
    % plot(x,mean(data.pei(:,:,izOut),2),'r--','Linewidth',1,'Marker','o');
    % hold on;
    % plot(x,mean(data.pei(:,:,izOut),2),'b--','Linewidth',1,'Marker','o');
    % title('P_e Profile');
    % xlabel('x');
    % ylabel('p_e');
    % xlim([0,alx]);
    % grid on;
    % set(gca,'xTick',x_sol);
    % set(gca,'LineWidth',1);
    % set(gca,'yGrid','off');

    % subplot(313)
    % plot(x,mean(data.vey(:,:,2),2),'r--','Linewidth',1,'Marker','o');
    % hold on;
    % plot(x,mean(data.vey(:,:,2),2),'b--','Linewidth',1,'Marker','o');
    % title('v_y Profile');
    % xlabel('x');
    % ylabel('v_y');
    % xlim([0,alx]);
    % grid on;
    % set(gca,'xTick',x_sol);
    % set(gca,'LineWidth',1);
    % set(gca,'yGrid','off');

    % filename='profile';
    % saveas(gcf,[filename,'.png']);
    % close;

    % plot phi profile
    figure;
    set(gca,'FontSize',14);
    plot(x,mean(data.phi(:,:,izOut),2),'r--','Linewidth',1,'Marker','o');
    hold on;
    plot(x,mean(data.phi(:,:,izIn),2),'b--','Linewidth',1,'Marker','o');
    title('Potential Profile');
    legend('Outside','Inside');
    xlabel('x');
    ylabel('phi');
    xlim([0,alx]);
    grid on;
    set(gca,'xTick',x_sol);
    set(gca,'LineWidth',1);
    set(gca,'yGrid','off');
    filename=sprintf('aprofile%4.4d',nt);
    % saveas(gcf,filename);     % for Matlab r2014a save as .fig
    saveas(gcf,[filename,'.png']);      % for Matlab r2014b save as .png
    close;

    % plot pei profile
    figure;
    set(gca,'FontSize',14);
    plot(x,mean(data.pei(:,:,izOut),2),'r--','Linewidth',1,'Marker','o');
    hold on;
    plot(x,mean(data.pei(:,:,izIn),2),'b--','Linewidth',1,'Marker','o');
    title('P_e Profile');
    legend('Outside','Inside');
    xlabel('x');
    ylabel('p_e');
    xlim([0,alx]);
    grid on;
    set(gca,'xTick',x_sol);
    set(gca,'LineWidth',1);
    set(gca,'yGrid','off');
    filename=sprintf('bprofile%4.4d',nt);
    % saveas(gcf,filename);     % for Matlab r2014a save as .fig
    saveas(gcf,[filename,'.png']);      % for Matlab r2014b save as .png
    close;

    % plot vey profile
    figure;
    set(gca,'FontSize',14);
    plot(x,mean(data.vey(:,:,izOut),2),'r--','Linewidth',1,'Marker','o');
    hold on;
    plot(x,mean(data.vey(:,:,izIn),2),'b--','Linewidth',1,'Marker','o');
    title('v_y Profile');
    legend('Outside','Inside');
    xlabel('x');
    ylabel('v_y');
    xlim([0,alx]);
    grid on;
    set(gca,'xTick',x_sol);
    set(gca,'LineWidth',1);
    set(gca,'yGrid','off');
    filename=sprintf('cprofile%4.4d',nt);
    % saveas(gcf,filename);       % for Matlab r2014a save as .fig
    saveas(gcf,[filename,'.png']);        % for Matlab r2014b save as .png
    close;
end

% % for Matlab r2014b
% % must load .fig and convert .png using normal for
% for nt=1:nts
%     filename=sprintf('%2.2dz%4.4d',izOut,nt);
%     fig=openfig(filename,'new','invisible');
%     saveas(fig,[filename,'.png']);
%     delete([filename,'.fig']);

%     filename=sprintf('%2.2dz%4.4d',izIn,nt);
%     fig=openfig(filename,'new','invisible');
%     saveas(fig,[filename,'.png']);
%     delete([filename,'.fig']);

%     filename=sprintf('%2.2dy%4.4d',iy,nt);
%     fig=openfig(filename,'new','invisible');
%     saveas(fig,[filename,'.png']);
%     delete([filename,'.fig']);

%     filename=sprintf('%2.2dx%4.4d',ix,nt);
%     fig=openfig(filename,'new','invisible');
%     saveas(fig,[filename,'.png']);
%     delete([filename,'.fig']);

%     filename=sprintf('aprofile%4.4d',nt);
%     fig=openfig(filename,'new','invisible');
%     saveas(fig,[filename,'.png']);
%     delete([filename,'.fig']);

%     filename=sprintf('bprofile%4.4d',nt);
%     fig=openfig(filename,'new','invisible');
%     saveas(fig,[filename,'.png']);
%     delete([filename,'.fig']);

%     filename=sprintf('cprofile%4.4d',nt);
%     fig=openfig(filename,'new','invisible');
%     saveas(fig,[filename,'.png']);
%     delete([filename,'.fig']);
% end

% delete(gcp('nocreate'));
close all