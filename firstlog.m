% by Mr. B. Li
function E=firstlog(A,dx)

dx=0.5/dx;

nx=length(A);
E=zeros(1,nx);
for i=2:nx-1
    E(i)= log(A(i+1)/A(i-1))*dx;
end
E(nx)= E(nx-1);
E(1)= E(2);

end