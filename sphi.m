%% boundary condition
pei=sbc(pei,2);
pei=sbc(pei,1,bc_p);%pi=0, free-fix
wi=sbc(wi,2);
deni=sbc(deni,2);
deni=sbc(deni,1,bc_n);
%% solve Phi
wm=zeros(nx0,ny0);%fix
phim=wm;
phi=zeros(size(phi));
%nx1=nx-1;cr=zeros(nx1,ny0);%free

for k=2:nz-1
    
    
    % cr(1:nx1,1:ny0)= -wi(1:nx1,2:ny0+1,k);%free-fix
    
    wm(1:nx0,1:ny0)= wi(2:nx0+1,2:ny0+1,k);% fix-fix
    wm(1,1:ny0)= wi(2,2:ny0+1,k); %p_left/dx^2;
    wm(nx0,1:ny0)= wi(nx0+1,2:ny0+1,k); %p_right/dx^2;
    wm=wm/rhoL^2;%180914 change normalization
    wm = fft(wm,[],2);
    
    for j=1:ny0
        % bi=-0.6*m2(j).^2./x.^4-0.7/dx^2.*(1./x1(2:end).^2+1./x1(1:end-1).^2);
        lap_op=poisson_A+diag(-0.6*m(j).^2./x(2:end-1).^4);%Laplace operator
        phim(:,j)=lap_op\wm(:,j);
    end
    
    phim = ifft(phim,[],2);
    
    phi(2:nx0+1,2:ny0+1,k)=real(phim(1:nx0,1:ny0));
    
end
phi=sbc(phi,2);
% debug(phi);%debug