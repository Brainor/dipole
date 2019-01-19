function sden(f,fi)
global deni den source_den cpe isextended isdeltaf den0
% \partial N/\partial t=curv(G)-convect(N)+S_N+dif(N)
if isdeltaf
    fxy=convect(deni,den0);
else
    fxy=convect(deni);
end
dif_den=dif(den);
fxy=fxy-dif_den-source_den;
if isextended;fxy=fxy-cpe;end
work=den*f+deni*fi-fxy;

if (f<0.6)
    den=deni;
    deni=work;
else
    deni=work;
end
end
