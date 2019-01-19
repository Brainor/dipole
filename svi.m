function svi(f,fi)
	global   vi vii  dpdz % nx ny nz pei pii petilde 

	%ptot=pei+pii;tpe=sbdel(ptot);
    %tpe=sbdel(petilde);
    %tpe=sbdel(pei);
	
	%fxy=convect(vii);
    dif2t=dif3(vi);
    %dif2t=dif4(vi);
	%fxy=fxy+tpe-dif2t;
    	%fxy=tpe-dif2t;
        	fxy=dpdz-dif2t;

	%work = vii;
	%work(2:nx-1,2:ny-1,2:nz-1)=vi(2:nx-1,2:ny-1,2:nz-1)*f+vii(2:nx-1,2:ny-1,2:nz-1)*fi-fxy(2:nx-1,2:ny-1,2:nz-1);

    work=vi*f+vii*fi-fxy;
    
	if(f<0.6)
		vi=vii;
		vii=work;
	else
		vii=work;
	end
end
