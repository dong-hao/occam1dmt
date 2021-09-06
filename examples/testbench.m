% simple testbench script, for 1D MT (magnetotellurics) occam inversion
% DONG Hao
% 2011/06/25
% Golmud
% ======================================================================= %
clear
addpath(genpath('..'),'-end');
% some settings here
% terminating RMS misfit
Trms=1.0;
% number of maximum iteration
Niter=20; 
% read a 24-layered model file
fid = fopen('aushield.mod');
tmp = textscan(fid,'%f %f','CommentStyle', '#');
fclose(fid);
res0 = tmp{2}';
l = tmp{1}(1:end-1)';
% depth of each layer INTERFACE, note that z1=0
z = cumsum([0 l]);
nz = length(l);
% read a MT sounding data file, which are the 
% 1D MT data published by Cull (1984)
fid = fopen('aushield.dat');
tmp = textscan(fid,'%f %f %f %f %f','CommentStyle', '#');
fclose(fid);
% the data are (literally) copied from Table 5 of the Constable 1987 paper
% see: 
% 
% Constable, S. C., Parker, R. L., & Constable, C. G. (1987). Occam’s 
% inversion: A practical algorithm for generating smooth models from 
% electromagnetic sounding data. Geophysics, 52(3), 289–300. 
% 
% Cull, J. P. (1985). Magnetotelluric soundings over a Precambrian contact 
% in Australia. Geophys. J. Roy. Astr. Sot., 80, 661-675.
period = tmp{1}';
rho0 = tmp{2}';
erho = tmp{3}'; 
phs0 = tmp{4}';
ephs = tmp{5}'; 

% a little setup to convert periods to frequency, and degree to rad 
freq = 1./period;
phs0 = phs0/180*pi;
ephs = ephs/180*pi;
[resi,rho,phs]=occam1dmt(z,res0,freq,rho0,erho,phs0,ephs,Trms,Niter);
figure(1);
l(end+1)=l(end)*1.5;
plotlayer_log(res0, l,'r');
hold on;
plotlayer_log(resi,l,'b');
figure(2);
plot1derr(freq,rho0,phs0,'rho','rx',erho,ephs);
plot1derr(freq,rho,phs,'rho','b-');
legend('obs','rsp');
% hasta la vista(?