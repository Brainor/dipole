function h=sp0(f)
% SP0 the turbulence. Minus the mean value in the y direction
h=zeros(size(f));
h(2:end-1,:,2:end-1)=f(2:end-1,:,2:end-1)-repmat(mean(f(2:end-1,:,2:end-1),2),[1,size(h,2),1]);
end