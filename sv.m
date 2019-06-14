function sv(f,fi)
    global vi vii
    dif2t=dif3(vi);
    fxy=dpdz-dif2t;
    
    work=vi*f+vii*fi-fxy;
    
    if(f<0.6)
        vi=vii;
        vii=work;
    else
        vii=work;
    end
end
