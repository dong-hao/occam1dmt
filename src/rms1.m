function rmsquare=rms1(obs,res,err)
% script to calculate Root Mean Square
% obs: observed
% res: responsed
% err: data *absolute* error
switch nargin
    case 0
        error('not enough input arguments, 2 at least')
    case 1
        error('not enough input arguments, 2 at least')       
    case 2
        err = ones(size(obs)); 
end
misfit=abs(obs-res)./err;
misfit=misfit.*misfit;
N=length(obs);
rmsquare=sqrt(sum(misfit)/N);
return