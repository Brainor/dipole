function h=curv(g)
%CURV curvature term
%%
% $$\frac{\partial\delta V}{\partial x}\frac{\partial p_e}{\partial y}=-4x^{4\gamma-5}\frac{\partial G}{\partial y}$$
global dyt2 x gamma
h=zeros(size(g));

h(2:end-1,2:end-1,2:end-1)=repmat(-4*x(2:end-1)'.^(4*gamma-5),[1,size(h,2)-2,size(h,3)-2]).*(g(2:end-1,3:end,2:end-1)-g(2:end-1,1:end-2,2:end-1))*dyt2;

