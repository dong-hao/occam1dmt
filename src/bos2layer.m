function [z0,res0]=bos2layer(rho0,phs0,freq)
% a silly function to guess a initial model from bostick transformation for
% 1D layered inversion
% first set up resistivity gauge
rmax=10;rmin=-10;
[resb,depth]=bostick(rho0,phs0,1./freq,'Bostick');
resb(end+1)=resb(end);
% try removing unpratical resistivity values first...
resb(resb>rmax)=rmax;
resb(resb<rmin)=rmin;
depth(end+1)=depth(end)*1.5;
% try to guess initial model from bostick transformation
% the first layer...
l1= ceil(depth(1)/10)*10;
le= ceil(depth(end)/100)*100*1.5;
ilayer=l1; n=1;
% maximum 100 layers
z0=zeros(1,100);res0=z0;
z0(1)=-1;res0(1)=0;
while z0(n)<=le
    n=n+1;
    ilayer=ilayer*1.2;
    z0(n)=z0(n-1)+ilayer;
    idx=istrap(depth,z0(n-1),z0(n));
    if isempty(idx)
        idx=find(depth>z0(n-1),1);
        if isempty(idx)
            idx=length(depth);
        else
            idx=find(depth>z0(n),1);
        end
    end
    res0(n)=mean(resb(idx));
end
idx=find(z0==0,1);
z0(1:idx-2)=z0(2:idx-1);
z0(idx-1:end)=[];
res0(1:idx-2)=res0(2:idx-1);
res0(idx-1:end)=[];
return 