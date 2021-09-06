function [res,rho,phs]=occam1dmt(z,res0,freq,rho0,erho,phs0,ephs,Trms,Niter)
% main function of occam inversion
% for inversion of log10 app. resistivity and impedance phase.
% please note APPARENT RESISTIVITY AND PHASE ARE INVERTED simutanously 
% for now. 
% DONG Hao
% 2011/06/25
% Golmud
% Note: this is but a toy scheme (no one is using 1D anymore, oh wait...) 
% when I was fiddling my 3D inversion code for my PhD. Thesis. 
% Although this sort of worked (as far as I recall), I have not extensively 
% tested this script - use with caution. 
%=========================================================================%
% input and output parametres:
% 
% z:        array of layer DEPTH of each layer INTERFACE of the starting
%           model, z1=0(surface of the earth)
% res0:     array of resistivity of each layer of the starting model
% freq:     array of input data frequency 
% rho0:     array of input data apparent resistivity
% erho:     array of input resistivity (absolute) error 
% phs0:     array of input data impedance phase
% ephs:     array of input phase (absolute) error 
% Trms:     target Root Mean Square misfit for iteration
% Niter:    maximum number of iteration
%=========================================================================%
% other parametres that might be useful:
%
% Nz:       number of model layers
% Nf:       number of frequencies 
% sigma:    conductivity, 1/res
% lambda0:  initial lagrange multiplier
%=========================================================================%
% setup some parametres
M=length(res0);
res=res0;
lambda=10;
dlambda=10;
%==================initialize inversion iteration=========================%
for iter=1:Niter
    % the outer loop
    disp(['================ iteration ', num2str(iter),' =================='])
%   plot_along(res,z,iter);
    lambdar=lambda*dlambda; % lambda on the right
    lambdal=lambda/dlambda; % lambda on the left
    [J1,J2,C]=jacob1(freq,res,z);
    rho=log10(7.8957E-06*abs(C(:,1)').^2.*freq);
    phs=atan2(real(C(:,1)),-imag(C(:,1)))';
    imrms=rms1([rho,phs],[rho0,phs0],[erho,ephs]);     % starting rms
    fprintf('! previous RMS = %5.3f \n', imrms);
    drho=rho0-rho;
    dphs=phs0-phs;
    resm=occam(drho,erho,dphs,ephs,res,lambda,J1,J2);  % resistivity in the middle
    resr=occam(drho,erho,dphs,ephs,res,lambdar,J1,J2); % resistivity on the right
    resl=occam(drho,erho,dphs,ephs,res,lambdal,J1,J2); % resistivity on the left
    [rhom,phsm] = mt1dfwd10(freq,resm,z,'rho');
    [rhor,phsr] = mt1dfwd10(freq,resr,z,'rho');
    [rhol,phsl] = mt1dfwd10(freq,resl,z,'rho');
    [ChiSl,ChiSm,ChiSr]=dispfit([rho0,phs0], [rhom,phsm], [rhol,phsl], ...
        [rhor,phsr],[erho,ephs], resl ,resm, resr, lambdal ,lambda, lambdar);  
    for ifind=1:10
        % the inner loop to sweep through different Lagrange multipliers
        % (or mu in origin Occam theory)
        if ChiSm <= ChiSl && ChiSm <= ChiSr
            disp('=============regional minimum found=============')
            break
        elseif ChiSm > ChiSl && ChiSr > ChiSl
            disp('================<< searching <<================')
            if lambdal/dlambda<=0.0001
                disp('! lambda too small, restart iteration')
                break
            end 
            resr=resm;
            lambdar=lambda;
            resm=resl;
            lambda=lambdal;
            lambdal=lambda/dlambda;
            resl=occam(drho,erho,dphs,ephs,res,lambdal,J1,J2);  % resistivity on the left
        elseif ChiSr < ChiSm && ChiSr < ChiSl
            disp('================>> searching >>================')
            if lambda*dlambda>=100000
                disp('! lambda too large, restart iteration')
                break
            end
            resl=resm;
            lambdal=lambda;
            resm=resr;
            lambda=lambdar;
            lambdar=lambda*dlambda;
            resr=occam(drho,erho,dphs,ephs,res,lambdar,J1,J2);  % resistivity on the right
        else 
            % seems not converging
            % use the middle value for next inversion iteration
            disp('=============cannot find a local minimum=============')
            break
        end
        [rhom,phsm] = mt1dfwd10(freq,resm,z,'rho');
        [rhol,phsl] = mt1dfwd10(freq,resl,z,'rho');
        [rhor,phsr] = mt1dfwd10(freq,resr,z,'rho');
        [ChiSl,ChiSm,ChiSr]=dispfit([rho0,phs0], [rhom,phsm], [rhol,phsl], ...
        [rhor,phsr],[erho,ephs], resl ,resm, resr, lambdal ,lambda, lambdar);  
    end
    res=resm';
    rho=rhom;
    phs=phsm;
    irms=rms1([rho, phs],[rho0, phs0],[erho, ephs]);
    if irms<=imrms
        % So far so good...
        disp(['! finishing iteration # ' num2str(iter)]);
    else
        % Huston, we have a convergence problem. 
        disp(['! warning: iteration # ' num2str(iter) ' not converged']);
        disp('! try stablizing inversion with smaller lambda...');
        disp(['! RMS= ' num2str(irms)]);
        lambda=lambda/dlambda.^2;
        continue
    end
    if irms<=Trms
        % check if the desired rms is reached
        fprintf('! current RMS = %5.3f \n', irms);
        fprintf('! target rms (%5.3f given by user) reached \n', Trms);
        disp('! try to find a smoothest model that fit the data')
        rmsr=irms;
        drho=rho0-rho;
        dphs=phs0-phs;
        while rmsr<=Trms
            disp('================>> searching >>================')
            lambda=lambda*sqrt(dlambda);
            resr=occam(drho,erho,dphs,ephs,res,lambda,J1,J2);
            [rhor,phsr] = mt1dfwd10(freq,resr,z,'rho');
            rmsr=rms1([rho0,phs0],[rhor,phsr],[erho,ephs]);
            fprintf('! evaluating a smoother model with Lambda = %5.3f \n', lambda);
            fprintf('! current RMS = %5.3f \n', rmsr);
            if rmsr < Trms 
                fprintf('! model accepted \n');
                rhom = rhor;
                phsm = phsr;
                resm = resr;
            else
                fprintf('! model rejected \n');
            end
        end
        disp('! picking up the smoothest model...')
        disp('! exiting...')
        res=resm';        
        break
    end
    if ifind>=10
        disp('! maximum lambda search reached, start next search...')
    end
end
% =========== end main iterations here =========== %
if iter>=Niter
    disp('! iteration limit reached, stop...')
else
    rho = rhom;
    phs = phsm;
end
frms=rms1([rho,phs],[rho0,phs0],[erho,ephs]);  
ChiS=chi2([rho,phs],[rho0,phs0],[erho,ephs]);
if ChiS<=2*M
    disp('! WARNING: desired chi2 reached. ');
    disp('! You might be over fitting the data! ')
end
fprintf('! final RMS = %5.3f \n', frms);
return

function res=occam(dpara1,epara1,dpara2,epara2,res,lambda,J1,J2)
% main function for a single (multi-variance) occam iteration 
% the dpara and epara be the diff and variance of the parameter, which can
% be apparent res, impedance phase, and (real or imag) part of complex
% impedance C.
N=size(J1,2);% number of model
D=mk_pmat(N); % roughness matrix
% penalty matrix
W1=diag(1./epara1);
W2=diag(1./epara2);
a=lambda*(D)'*D;%+diag(ones(N,1));
% construct weighted Jacobian premultiplied by TRANS(W.J)
wjtwj=(W1*J1)'*(W1*J1)+(W2*J2)'*(W2*J2);
% construct weighted translated data premultiplied by TRANS(W.J)
wjtwd=(W1*J1)'*(W1*(dpara1'+J1*res'))+(W2*J2)'*(W2*(dpara2'+J2*res'));
res=(a+wjtwj)\wjtwd;
return

function DEL=mk_pmat(N)
% constructing roughness matrix 
% for L1 norm
DEL=diag(ones(N,1),0)-diag(ones(N-1,1),-1); 
DEL(1,1)=0;
return 

function [ChiSl,ChiSm,ChiSr]=dispfit(co, cm, cl, cr, stderr, resl, resm, resr, lambdal ,lambda, lambdar)
ChiSm=chi2(co,cm,stderr);
ChiSl=chi2(co,cl,stderr);
ChiSr=chi2(co,cr,stderr);
roughl=roughness1(resl,1);
roughm=roughness1(resm,1);
roughr=roughness1(resr,1);
fprintf('               L          M          R\n');
fprintf('Chi Square   : %8.4e %8.4e %8.4e \n',...
    ChiSl, ChiSm, ChiSr);
fprintf('Roughness(L1): %8.4e %8.4e %8.4e \n',...
    roughl, roughm, roughr);
fprintf('Lambda       : %8.5g %8.5g %8.5g \n',...
    lambdal, lambda, lambdar );
return