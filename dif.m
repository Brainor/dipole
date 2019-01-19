function h=dif(f)
%   $$D\nabla_\perp^2 f$$
	global  dxt2d dyt2d
	h=zeros(size(f));
    h(2:end-1,2:end-1,2:end-1) =(f(3:end,2:end-1,2:end-1)-2*f(2:end-1,2:end-1,2:end-1)+f(1:end-2,2:end-1,2:end-1))*dxt2d...
                              +(f(2:end-1,3:end,2:end-1)-2*f(2:end-1,2:end-1,2:end-1)+f(2:end-1,1:end-2,2:end-1))*dyt2d;
end
