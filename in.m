%%%%%%%%%%
% Update Log
% 180830 Coefficient of Poisson Equation is wrong
% 180914 Change w normalization
% 181226 New equations
% 
%%%%%%%%%%
clear variables global; close all;
addpath(genpath([regexp(pwd,'^.*data','match','once'),'/../code/190119.deltaf.owk']));%use '/' to indicate seperator of path, compatable with Windows and UNIX
if ~isfolder('data');mkdir('data');end%store data
if ~isfolder('plot');mkdir('plot');end%store plot

global nx0 ny0 nz  dif tau gamma rhoL bc_p bc_n isextended isdeltaf
restart=0;
%restart=1;
isextended=1;%ideal 0; extended 1
isdeltaf=1;%using deltaf
bc_p=[0,0];bc_n=[0,0];
rhoL=.04;%rho_s/L

% time grid.
% linear Loop ntp times to produce 1 plot. Loop nts*ntp times in total.
ntp=100;  nts=100;
% box grid number
nx0=64;nx0=50*3;
ny0=160;ny0=128*3;
% time step
tau=.002;%linear
tau=.005;
tau=.0005/2;
% difusion parameter
dif=0.01;
% source amplitude
sn=0.14;%linear,close
sp=0.14;
sn=0;sp=0;%no source
% initial profile
amp=1;
pert=1.e-5;
xs=1;xw=0.1;
%equilibrim term



% fixed parameters
gamma=5/3;
L0=0.9;
xmin=L0/2.5; xmax=L0/0.6;
alx=xmax-xmin;
aly=2*pi;
nz=3;%2d
% nz=18; %3d

main;
show_evo;% 2D profile (plot and video) in every ntp times. Profile at last slide at y=pi
show_profile_meanty;% profile with average of y and t
show_profile_meany; % profile with average of y
show_phase_flux_phi;% phase analysis of p_e and Phi
show_spectrum; % spectrum analysis


