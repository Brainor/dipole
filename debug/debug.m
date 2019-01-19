function debug(x)
%DEBUG show the max and min position and value of variable x
%   x is a 3-dimention-matrix
x=abs(x);
[M1,f]=min(x(:));
[f1{1},f1{2},f1{3}]=ind2sub(size(x),f);
[M2,f]=max(x(:));
[f2{1},f2{2},f2{3}]=ind2sub(size(x),f);
fprintf('%s:[(%d,%d,%d):%.2g;(%d,%d,%d):%.2g].',inputname(1),f1{:},M1,f2{:},M2);
end

