function plot1derr(freq,var1,var2,opt,ccode,err1,err2)
% a barn door function to plot 1d apparent resistivities and impedance
% phases...
% 
switch nargin
    case 0
        error('not enough input arguments, 3 at least')
    case 1
        error('not enough input arguments, 3 at least')       
    case 2 
        error('not enough input arguments, 3 at least')   
    case 3
        ccode='b-';
        opt = 'rho';
    case 4
        ccode='b-';
    case 6
        error('not enough input arguments, need 2 vars and 2 errs')  
end
mju0=4.0*pi*1E-7; 
omega=2*pi*freq; 
if strcmp(opt,'imped')==1
    rho=(var1.^2+var2.^2)*mju0.*omega';
    rho = log(rho);
    phase=atan2(var2,var1);
    rhoe=2*sqrt(err1.^2+err2.^2)./((var1.^2+var2.^2).^0.5);
    rhoe = 0.434*rhoe;
    phse=sqrt(err1.^2+err2.^2)./((var1.^2+var2.^2).^0.5);
else % in log10 domain
    rho=var1;
    phase=var2;
    if nargin > 5
        rhoe=err1;
        phse=err2;
    end
end
a1=subplot(2,1,1);
if nargin <= 5
    plot(a1,freq,rho,ccode);
else
    errorbar(a1,freq,rho,rhoe,ccode);
end
% loglog(a1,freq,rho,mark);
hold(a1,'on');
set(a1,'xdir','reverse')
set(a1,'xscale','log');
set(a1,'xgrid','on','ygrid','on')
a2=subplot(2,1,2);
phase=phase/pi*180;
if nargin <= 5
    plot(a2,freq,phase,ccode);
else
    phse = phse/pi*180;
    errorbar(a2,freq,phase,phse,ccode);
end
hold(a2,'on');
set(a2,'xdir','reverse','ylim',[0 90],'ytick',[0 15 30 45 60 75 90]);
set(a2,'xscale','log');
set(a2,'xgrid','on','ygrid','on')
return