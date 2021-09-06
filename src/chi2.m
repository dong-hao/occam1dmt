function X2=chi2(d0,d,stderr)
% script to calculate Chi square statistic
%=========================================================================%
% obs: observed
% res: responsed
% stderr: standard variance of the data
% c2: chi square (normalized)
% expected Chi square should be 2*length(obs)
if size(d0)~=size(d) 
    error('size of the input vectors must be identical')
elseif size(d0)~=size(stderr)
    error('size of the input vectors must be identical')
end
W=diag(1./stderr);
X2=sum((W*(d0'-d')).^2);
return
