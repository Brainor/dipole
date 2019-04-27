global xmax xmin x x1
%参数设置
isextended=0;%ideal:0, extended:1
m=4;
%默认参数
gamma=5/3;
L0=0.9;
xmin=L0/2.5; xmax=L0/0.6;
nx0=1000;
rhoL=.5;%rho_s/L

alx=xmax-xmin;
nx=nx0+2;%nx positions in total including the end-points, nx0 not including
dx=alx/(nx-1);
x=xmin:dx:xmax;
omegad=4./x(2:end-1);

%平衡剖面设置

den0=h_profile(0.1,0.5);%den0=ones(size(x));
pe0=h_profile(0,1);

den1=den0(2:end)-den0(1:end-1);%hN(0.5):dx:hN(N+0.5)
pe1=pe0(2:end)-pe0(1:end-1);%hG(0.5):dx:hG(N+0.5)
x1=xmin+dx/2:dx:xmax-dx/2;%x(0.5):dx:x(N+0.5)
den1=h_profile(0.1,0.5,1);%den1=ones(size(x1));
pe1=h_profile(0,1,1);

%先记录ideal
M=zeros(2*nx0);
N=zeros(2*nx0);

M(1:nx0,1:nx0)=eye(nx0);

ai=0.7/dx^2.*x1(2:end-1).^-2.*den1(2:end-1);
bi=-0.7/dx^2.*(x1(2:end).^-2.*den1(2:end)+x1(1:end-1).^-2.*den1(1:end-1))-0.6*m^2*den0(2:end-1).*x(2:end-1).^-4;
ci=0.7/dx^2.*x1(2:end-1).^-2.*den1(2:end-1);
M(nx0+1:2*nx0,nx0+1:2*nx0)=diag(bi)+diag(ai,-1)+diag(ci,1);

N(1:nx0,nx0+1:2*nx0)=m/dx*diag((pe1(2:end)-pe1(1:end-1)).*x(2:end-1).^(4*gamma));
% N(nx0+1:2*nx0,1:nx0)=m/rhoL^2*diag(omegad.*x(2:end-1).^(-4*gamma));
N(nx0+1:2*nx0,1:nx0)=m/rhoL^2*diag(4*x(2:end-1).^(-5));

if isextended
    M(nx0+1:3*nx0,nx0+1:3*nx0)=M;
    N(nx0+1:3*nx0,nx0+1:3*nx0)=N;
    N(1:nx0,nx0+1:2*nx0)=m*omegad.*diag(x(2:end-1).^(-4));
    N(1:nx0,2*nx0+1:3*nx0)=m/dx*diag(den1(2:end)-den1(1:end-1));
    N(nx0+1:2*nx0,1:nx0)=-m*gamma*diag(omegad.*x(2:end-1).^(8*gamma-4).*(pe0(2:end-1)./den0(2:end-1)).^2);
    N(nx0+1:2*nx0,nx0+1:2*nx0)=2*m*gamma.*diag(omegad.*x(2:end-1).^(4*gamma-4).*(pe0(2:end-1)./den0(2:end-1)));
end
M=sparse(M);
N=sparse(N);

[V,w]=eigs(N,M,3,'largestimag');
%%
draw_plot({L0./x(2:end-1),real(V(nx0+1:2*nx0,:))},['$$m=',num2str(m),'$$, $$\gamma=',num2str(imag(w(1))),'$$'],'L','');


function h=h_profile(h0,c,~)
global xmax xmin x x1
b=c*(xmax-1)/(1-xmin);
a=(1-h0)/(xmax-1)^b/(1-xmin)^c;
if nargin==2
    h=h0+a*(xmax-x).^b.*(x-xmin).^c;
else
    h=h0+a*(xmax-x1).^b.*(x1-xmin).^c;
end
end

