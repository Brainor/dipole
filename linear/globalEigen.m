global xmax xmin x_eigen x1_eigen
%参数设置
isextended=0;%ideal:0, extended:1
isdiffusive=1;%考虑耗散项
D=.06;
ms=1:10;growthrate=zeros(1,length(ms));phim_eigen=[];
for i=1:length(ms)
    m=4;
    m=ms(i);
    %默认参数
    gamma=5/3;
    L0=0.9;
    xmin=L0/2.5; xmax=L0/0.6;
    nx0=500*2*2*2;
    rhoL=2;%rho_s/L
    rhoL=.04;
    
    alx=xmax-xmin;
    dx_eigen=alx/(nx0+1);
    x_eigen=xmin:dx_eigen:xmax;
    omegad=4./x_eigen(2:end-1);
    
    %平衡剖面设置
    
    den0=h_profile(0.1,0.5);%den0=ones(size(x));
    pe0=h_profile(0,1);
    
    den1=den0(2:end)-den0(1:end-1);%hN(0.5):dx:hN(N+0.5)
    pe1=pe0(2:end)-pe0(1:end-1);%hG(0.5):dx:hG(N+0.5)
    x1_eigen=xmin+dx_eigen/2:dx_eigen:xmax-dx_eigen/2;%x(0.5):dx:x(N+0.5)
    den1=h_profile(0.1,0.5,1);%den1=ones(size(x1));
    pe1=h_profile(0,1,1);
    
    % den0=ones(size(den0));
    % den1=ones(size(den1));
    % %ivp形状
    % den0=amp*exp(-((x-xs_n)/xw).^2);
    % pe0=amp*exp(-((x-xs_p)/xw).^2);
    % den1=amp*exp(-((x1-xs_n)/xw).^2);
    % pe1=amp*exp(-((x1-xs_p)/xw).^2);
    
    %先记录ideal
    M=zeros(2*nx0);
    N=zeros(2*nx0);
    
    M(1:nx0,1:nx0)=eye(nx0);
    
    ai=0.7/dx_eigen^2.*x1_eigen(2:end-1).^-2.*den1(2:end-1);
    bi=-0.7/dx_eigen^2.*(x1_eigen(2:end).^-2.*den1(2:end)+x1_eigen(1:end-1).^-2.*den1(1:end-1))-0.6*m^2*den0(2:end-1).*x_eigen(2:end-1).^-4;
    ci=0.7/dx_eigen^2.*x1_eigen(2:end-1).^-2.*den1(2:end-1);
    M(nx0+1:2*nx0,nx0+1:2*nx0)=diag(bi)+diag(ai,-1)+diag(ci,1);
    
    N(1:nx0,nx0+1:2*nx0)=m/dx_eigen*diag((pe1(2:end)-pe1(1:end-1)).*x_eigen(2:end-1).^(4*gamma));
    % N(nx0+1:2*nx0,1:nx0)=m/rhoL^2*diag(omegad.*x(2:end-1).^(-4*gamma));
    N(nx0+1:2*nx0,1:nx0)=m/rhoL^2*diag(omegad.*x_eigen(2:end-1).^(-4));
    
    if isextended
        M(nx0+1:3*nx0,nx0+1:3*nx0)=M;
        N(nx0+1:3*nx0,nx0+1:3*nx0)=N;
        N(1:nx0,nx0+1:2*nx0)=m*omegad.*diag(x_eigen(2:end-1).^(-4));
        N(1:nx0,2*nx0+1:3*nx0)=m/dx_eigen*diag(den1(2:end)-den1(1:end-1));
        N(nx0+1:2*nx0,1:nx0)=-m*gamma*diag(omegad.*x_eigen(2:end-1).^(8*gamma-4).*(pe0(2:end-1)./den0(2:end-1)).^2);
        N(nx0+1:2*nx0,nx0+1:2*nx0)=2*m*gamma.*diag(omegad.*x_eigen(2:end-1).^(4*gamma-4).*(pe0(2:end-1)./den0(2:end-1)));
    end
    M=sparse(M);
    N=sparse(N);
    
    [V,w]=eigs(N,M,3,'largestimag');
    if isdiffusive
        growthrate(i)=imag(w(1))-D*m^2;
    else
        growthrate(i)=imag(w(1));
    end
    phim_eigen(i,:)=real(V(nx0+1:end,1));
end
draw_plot({ms,growthrate},{'growth rate',sprintf('$$D=%g$$',D)},'m','\gamma');
%%
subplot(1,2,1)
% draw_plot({L0./x(2:end-1),real(V(nx0+1:2*nx0,:))},['$$m=',num2str(m),'$$, $$\rho_\star=',num2str(rhoL),'$$'],'L','\Re\{\Phi\}');
draw_plot({1./x_eigen(2:end-1),real(V(nx0+1:2*nx0,:))},['$$m=',num2str(m),'$$, $$\rho_\star=',num2str(rhoL),'$$'],'x','\Re\{\Phi\}','XAxisLocation','origin');
lgd=legend(cellfun(@num2str,num2cell(imag(diag(w))),'UniformOutput',false),'Location','best','Box','off');
lgd.Title.String='$$\gamma$$';lgd.Title.Interpreter='latex';
subplot(1,2,2)
% draw_plot({L0./x(2:end-1),[den0(2:end-1);pe0(2:end-1)]},'Equilibruim Profile','L');legend('$$h_N$$','$$h_G$$','Interpreter','latex','Box','off');
draw_plot({x_eigen(2:end-1),[den0(2:end-1);pe0(2:end-1)]},'Equilibruim Profile','x');legend('$$h_N$$','$$h_G$$','Interpreter','latex','Box','off');
% h=subplot(2,2,3);
X=x_eigen(2:end-1)'*cos(linspace(0,2*pi,100));
Y=x_eigen(2:end-1)'*sin(linspace(0,2*pi,100));
phi=real(V(nx0+1:2*nx0,1));
C=phi*cos(m*linspace(0,2*pi,100));
% draw_pcolor(h,{C,X,Y});
% streamline(X,Y,


%% Mike paper



function h=h_profile(h0,c,~)
global xmax xmin x_eigen x1_eigen
b=c*(xmax-1)/(1-xmin);
a=(1-h0)/(xmax-1)^b/(1-xmin)^c;
if nargin==2
    h=h0+a*(xmax-x_eigen).^b.*(x_eigen-xmin).^c;
else
    h=h0+a*(xmax-x1_eigen).^b.*(x1_eigen-xmin).^c;
end
end

