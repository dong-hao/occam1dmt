function [J1,J2,C]=jacob1(freq,res,z)
% rewritten Jacobian calculation script for 1D MT 
% for drhoa/dlog10(res) and dphs/dlog10(res)
% used for occam 1D MT inversion
% DONG Hao
% 2011/06/25
% Golmud
res=10.^res;
mju0=4.0*pi*1E-7; % Vacuum Permeability, unit is Vs/(Am) in SI.
sigma=1./res;     % conductivity
omega=2*pi*freq;  % angular frequency
NL=length(z);     % number of layers
NF=length(freq);  % number of freqs
J1=ones(NF, NL);   % jacobian for df/dx at each layer and each frequency.
J2=J1;
l=diff(z);
C=zeros(NF, NL);  
dCdr=ones(NF, NL);
q=zeros(NF,1);
for ifreq=1:NF % looping through frequencies
    % calculating the last layer first
    IOM=1j*omega(ifreq)*mju0;
    q(NL)=sqrt(IOM*sigma(NL)); 
    C(ifreq,NL)=1/q(NL);
    dCdr(ifreq,NL)=1/(q(NL)/sigma(NL)*2); % Note this is right only at the last layer
    if NL==1 % in the case of a homogenous half space
        continue
    end
    for ilayer=NL-1:-1:1
        q(ilayer)=sqrt(IOM*sigma(ilayer));
        qp1 = C(ifreq,ilayer+1)*q(ilayer)+1;
        qm1 = C(ifreq,ilayer+1)*q(ilayer)-1;
        theexp=exp(-2.*q(ilayer)*l(ilayer));
        opexp = theexp*qm1/qp1;
        omexp = 1 - opexp;
        opexp = 1 + opexp;
        dCdr(ifreq,ilayer) = 2*theexp*(C(ifreq,ilayer+1)/(qp1.^2)-(l(ilayer)*qm1/qp1))/omexp.^2;
        dCdr(ifreq,ilayer) = opexp/omexp/q(ilayer)/2 - dCdr(ifreq,ilayer);
        dCdr(ifreq,ilayer) = dCdr(ifreq,ilayer)*sigma(ilayer);
        % calculating dCi/dCi+1
        dciip1 = 4*theexp/(qp1*omexp).^2;
        % now calculating C using Wait's Recursion Formula  
        C(ifreq,ilayer)=opexp/omexp/q(ilayer);
        % cumproduction for dCi/dCi+1
        for klayer=ilayer+1:NL
            % note the J(dCdr) here is "drhoa/dlog10(res)"
            dCdr(ifreq,klayer)=dCdr(ifreq,klayer)*dciip1;
        end
    end
end
for ifreq=1:NF
    for ilayer=1:NL
% for log10(res) (logrithum)       
         J1(ifreq,ilayer)=2*res(ilayer)*(real(C(ifreq,1))*real(dCdr(ifreq,ilayer))+...
             imag(C(ifreq,1))*imag(dCdr(ifreq,ilayer)))/abs(C(ifreq,1)).^2;
% for res (linear)
%          J1(ifreq,ilayer)=omega(ifreq)*mju0*(real(C(ifreq,1))*real(dCdr(ifreq,ilayer))+...
%              imag(C(ifreq,1))*imag(dCdr(ifreq,ilayer)));
% for phs (linear in RAD)
         J2(ifreq,ilayer)=(real(C(ifreq,1))*imag(dCdr(ifreq,ilayer))-...
             imag(C(ifreq,1))*real(dCdr(ifreq,ilayer)))/abs(C(ifreq,1)).^2;
    end
end
return
