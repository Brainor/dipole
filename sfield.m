global cpe
cpe=curv(pei);

vey(2:end-1,2:end-1,2:end-1)=-dx2*(phi(3:end,2:end-1,2:end-1)-phi(1:end-2,2:end-1,2:end-1));%-\partial\Phi/\partial x
vey=sbc(vey,2);

%x_right is fix, set to 0;,x_left is free, set from differential equation
vex(2:end-1,2:end-1,2:end-1)=dy2*(phi(2:end-1,3:end,2:end-1)-phi(2:end-1,1:end-2,2:end-1)); %\partial\Phi/\partial y
% vex=sbcy(vex);%no need of BC?


