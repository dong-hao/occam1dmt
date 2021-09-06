function [o1,o2] = mt1dfwd10(freq,res,z,opt)
% a barn door 1d layered forward routine for 
% magnetotelluric 1D inversion
% log10 version (both the input resisitivity and output app. res will be in
% log10 space.
% apparent resistivity and phase are generated from layerd model.
% DONG Hao
% 2011/06/25
% Golmud
%=========================================================================%
% description of parameters:
% 
% freq:     array of frequency for response to be generated 
% sigma:    conductivity 1/res
% res:      array of resistivity of each layer
% z:        array of layer DEPTH of each layer INTERFACE, z1=0(surface of 
%           the earth)
% mju0:     absolute magnetic permeability of vacuum
% omega:    Angular frequency (2 pi frequency)
% opt:      option of outputing the forward, could be "rho" or "imped"
%=========================================================================%
% version 0.2
% DONG Hao
%=========================checking parameters=============================%
switch nargin
    case 0
        error('not enough input arguments, 2 at least')
    case 1
        error('not enough input arguments, 2 at least')       
    case 2
        z = 0; % treat the earth as half space
    case 3
        opt = 'imped';
end
if (size(z)~=size(res))
    if length(z)==length(res)
        res=res';
    else
	    disp('please check the input parametres ')
        error('res and z should have same size');
    end
end
%=========================================================================%
mju0=4.0*pi*1E-7; % unit is Wb/(A��m) in SI. 
omega=2*pi*freq; 
sigma=1./(10.^res);
NL=length(z); % number of layers
NF=length(freq); % number of freqs
C=zeros(NL, NF);% transfer function for each layer and each frequency.
q=zeros(1,NL);
for ifreq=1:NF
    q(NL)=sqrt(1j*omega(ifreq)*mju0*sigma(NL)); % calculating the last layer first
    C(NL,ifreq)=1/q(NL); % Note this is right only at the last layer
    if NL==1 % in the case of a homogenous half space
        continue
    end
    for ilayer=NL-1:-1:1
        q(ilayer)=sqrt(1j*omega(ifreq)*mju0*sigma(ilayer));
        l=z(ilayer+1)-z(ilayer); % depth of the layer
        % using Wait's Recursion Formula
        C(ilayer,ifreq)=1/q(ilayer)*(q(ilayer)*C(ilayer+1,ifreq)+tanh(q(ilayer)*l))...
            /(1+q(ilayer)*C(ilayer+1,ifreq)*tanh(q(ilayer)*l));
    end
end 
rho=C(1,:).*conj(C(1,:)).*omega*mju0; % only take the transfer function from the first layer
phase=atan2(real(C(1,:)),-imag(C(1,:)));
switch opt
    case 'rho'
        o1=log10(rho);
        o2=phase;
    case 'imped'
        Z=C(1,:).*(1j*omega);
        o1=real(Z);
        o2=-imag(Z);
end
return