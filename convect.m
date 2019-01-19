function h=convect(g,varargin)
%CONVECT convection term
% $$\partial v_ex g/\partial x+\partial v_ey g/\partial y=
% \partial/\partial x (g \partial\Phi/\partial y)-
% \partial/\partial y (g \partial\Phi/\partial x)
global  vex vey dyt2 dxt2
h = zeros(size(g));
h(2:end-1,2:end-1,2:end-1) = (vex(3:end,2:end-1,2:end-1).*g(3:end,2:end-1,2:end-1)-vex(1:end-2,2:end-1,2:end-1).*g(1:end-2,2:end-1,2:end-1))*dxt2...
    +(vey(2:end-1,3:end,2:end-1).*g(2:end-1,3:end,2:end-1)-vey(2:end-1,1:end-2,2:end-1).*g(2:end-1,1:end-2,2:end-1))*dyt2;
if(nargin==2)
    %using deltaf solution, with the parameters being g\delta,
    %add \partial g0/\partial x \partial \Phi/\partial y term
    %g0 is the equilibrium term of delta g
    g0=varargin{1};
    h(2:end-1,2:end-1,2:end-1)=h(2:end-1,2:end-1,2:end-1)+vex(2:end-1,2:end-1,2:end-1).*(g0(3:end,2:end-1,2:end-1)-g0(1:end-2,2:end-1,2:end-1))*dxt2;
end
