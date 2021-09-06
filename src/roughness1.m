function r=roughness1(sigma,opt)
% description of parameters:
%
% sigma: conductivity 1/res
% opt: option of outputing the forward, could be 1 for L1 and 2 for L2
%===============================================================%
% version 0.2
% DONG Hao
%====================checking parameters========================%
switch nargin
    case 0
        error('not enough input arguments, 1 at least')
    case 1
        opt = 1; % use L1 roughness
end
n=length(sigma); %number of layers 
switch opt
    case 1 %gradient
        r=abs(sigma(1:n-1)-sigma(2:n));
        r=r'*r;
    case 2 %laplacian
        r=abs(-sigma(1:n-2)+2*sigma(2:n-1)-1*sigma(3:n));
        r=r'*r;
end
return