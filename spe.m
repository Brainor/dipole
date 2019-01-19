function spe(f,fi)
global   pe pei source_p deni gamma isextended isdeltaf pe0
% \partial G/\partial t=\gamma curv(G^2/N)-convect(G)+S_G+dif(G)

if(isdeltaf)
    fxy=convect(pei,pe0);
else
    fxy=convect(pei);
end
dif_pe=dif(pe);
fxy=fxy-dif_pe-source_p;
if isextended;fxy=fxy-gamma*curv(2*pei-deni);end%¼ò»¯
work=pe*f+pei*fi-fxy;
% debug(pei);debug(deni);debug(pei.^2./deni);fprintf(newline);%debug
if(f<0.6)
    pe=pei;
    pei=work;
else
    pei=work;
end

end
