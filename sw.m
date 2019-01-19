function sw(f,fi)
global wi w cpe
% $\partial w/\partial t=curv(G)-connvect(w)+dif(w)$
fxy=convect(wi);
difw=dif(w);
fxy=-cpe+fxy-difw;
work=w*f+wi*fi-fxy;
if(f<0.6)
    w=wi;
    wi=work;
else
    wi=work;
end
end
