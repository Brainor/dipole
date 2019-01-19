% calc confinement time
% to de sent

na=1;
nb=10;
petotal=[];

for nt=na:nb
    load(['dat',sprintf('%4.4d',nt)]);
    % because of the B.C., choose x 1:end, y 1:end-1, z 2:end
    pet=pei(2:nx,1:end-1,2:end);
    
    
    % mean in y and z direction,and integrate in x direction
    petotal(end+1)=sum(mean(mean(pet,3),2));%*dx
end
% mean in time
peAVG=mean(petotal);


nz=2;%2d
   peAVG= sum(mean(pei(2:nx-1,:,nz),2));
    
SpAVG=sum(source_p(2:nx-1))/tau; %*dx

tau_cnfn=peAVG/SpAVG

% str1=['the confinement time is ',num2str(tau_cnfn)];
% msgbox(str1,'result');