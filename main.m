global  nx ny nz wi w vex vey  phi den deni pe pei ...
    dy2 dx2 dyt2 dxt2 dxt2d dyt2d source_p source_den ...
    x dtz m
%% Initial Value Setting
dtz=0;
nx=nx0+2;%nx positions in total including the end-points, nx0 not including
dx=alx/(nx-1);
x=xmin:dx:xmax;
ny=ny0+2;%ny0 positions in a circle, ny adds first and end points as boundary condition
dy=aly/ny0;%periodical condition, so it's a loop

% differential operators
dxt2=tau/(2.*dx);%dt/2dx
dyt2=tau/(2.*dy);
dx2=1./(2.*dx);%1/2dx
dy2=1./(2.*dy);
dxt2d=tau/(dx^2)*dif;%dt/2dx*dif
dyt2d=tau/(dy^2)*dif;

% source
s_p=tau*sp*exp(-((x-xs_p)/xw).^2);%pe source
s_den=tau*sn*exp(-((x-xs_n)/xw).^2);%n source
source_p = repmat(s_p',[1,ny,nz]);
source_den = repmat(s_den',[1,ny,nz]);

% equilibrim term
if(isempty(den0));den0=bc_n(1)+amp*exp(-((x'-xs)/xw).^2);end
if(isempty(pe0));pe0=amp*exp(-((x'-xs)/xw).^2);end
den0([1,end])=bc_n;den0=repmat(den0',[1,ny,nz]);%平衡剖面为沿x方向, 扩充成三维矩阵
den0=sbc(den0,2);den0=sbc(den0,3);
pe0([1,end])=bc_p;pe0=repmat(pe0',[1,ny,nz]);pe0=sbc(pe0,2);pe0=sbc(pe0,3);

% Poisson equation's matrix setting
m=[0:ny0/2,floor(-ny0/2)+1:-1];%m is the mode number in y direction. Phi is real, so m is symmetric at ny0/2. e.g. [0,1,2,-2,-1], [0,1,2,-1]
x1=xmin+dx/2:dx:xmax-dx/2;%x(0.5):dx:x(N+0.5)
ai=0.7/dx^2./x1(2:end-1).^2;
bi=-0.7/dx^2.*(1./x1(2:end).^2+1./x1(1:end-1).^2);
ci=0.7/dx^2./x1(2:end-1).^2;
poisson_A=diag(bi)+diag(ai,-1)+diag(ci,1);
poisson_A=sparse(poisson_A);

% Initialization
wi=zeros(nx,ny,nz);
phi=wi;vex=wi;vey=wi;pei=wi;deni=wi;
switch restart
    case 0
        pei(2:end-1,2:end-1,2:end-1)=pert*rand(nx-2,ny-2,nz-2);%deltape initial profile
        deni(2:end-1,2:end-1,2:end-1)=pert*rand(nx-2,ny-2,nz-2);%deltan initial profile
        if (~isdeltaf)
            pei(2:end-1,2:end-1,2:end-1)=pei(2:end-1,2:end-1,2:end-1)+pe0(2:end-1,2:end-1,2:end-1);
            deni(2:end-1,2:end-1,2:end-1)=deni(2:end-1,2:end-1,2:end-1)+den0(2:end-1,2:end-1,2:end-1);
            pei=sbc(pei,1,bc_p);deni=sbc(deni,1,bc_n);
        end
        pei=sbc(pei,3);pei=sbc(pei,2);
        deni=sbc(deni,3);deni=sbc(deni,2);
        den=deni;pe=pei;w=wi;
    case 1
        load rest.mat
end

draw_plot({squeeze(pei(:,floor(ny/2),2)),'-o'},'initial $$p_e$$ profile at $$y=\pi$$','x','p_e');
print('plot/prof_pe_initial','-dpng')
close

%% Solve Equations
for nt=1:nts
    fprintf('nt=%d/%d\n',nt,nts);
    for ntt=1:ntp
        %         fprintf('nt=%d/%d\n',ntt,ntp);%debug
        %         debug(wi);debug(pei);debug(deni);debug(w);debug(pe);debug(den);debug(phi);fprintf(newline);
        sfield
        f=0.5;
        fi=0.5;
        %svi(f,fi)
        sw(f,fi)
        spe(f,fi)
        sden(f,fi)
        sphi
        %         fprintf('halfway\n');%debug
        sfield
        f=1.0;
        fi=0.0;
        %svi(f,fi)
        sw(f,fi);
        spe(f,fi)
        sden(f,fi)
        sphi
        
        done=isfinite(pei(2,2,2));
        if done==0
            error('nan')
        end
    end
    save(sprintf('data/dat%4.4d',nt),'wi','pei','deni','phi','w','pe','den', 'vey','vex')
end

