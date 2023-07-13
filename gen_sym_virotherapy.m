%% Virotherapy model
% generate right-hand side and derivatives from symbolic expressions
%
% From: Chapter 4 Bruno thesis
%
%% Differential equations
%
% $$\frac{dx}{dt}=rx\left(1-\frac{x+y}{K}\right)-\beta xv,$$
% $$\frac{dy}{dt}=\beta x(t-\tau)v(t-\tau)-\mu yz-\delta y,$$
% $$\frac{dv}{dt}=b\delta y-\beta xv-\theta vz-\gamma v,$$
% $$\frac{dz}{dt}=szy-cz^2.$$
%
% Parameters are (in this order) |r|, |K|, |beta|, |mu|, |delta|, |b|, |theta|, |gamma|, |s|, |c|, |tau|
%% Add paths and load sym package if GNU Octave is used
clear
ddebiftoolpath='../../';
addpath(strcat(ddebiftoolpath,'ddebiftool'),...
strcat(ddebiftoolpath,'ddebiftool_extra_symbolic'));
if dde_isoctave()
pkg load symbolic
end
%% Create parameter names as strings and define fixed parameters
% The demo has the parameters |r|, |K|, |beta|, |mu|, |delta|, |b|, |theta|, |gamma|, |s|, |c| and |tau|
parnames={'r','K','beta','mu','delta','b','theta','gamma','s','c','tau'};
%% Create symbols for parameters, states and delayed states
% The array |par| is the array of symbols in the same order as parnames.
% Due to the following two lines we may, for example, use either tau or
% par(11) to refer to the delay.
syms(parnames{:}); % create symbols for the parameters.
par=cell2sym(parnames); % now r is par(1) etc
%% Define system using symbolic algebra
syms x xtau y ytau v vtau z ztau % create symbols for x(t) x(t-tau), y(t), y(t-tau), v(t), v(t-tau), z(t), z(t-tau)
dx_dt=r*x*(1-(x+y)/K)-beta*x*v;
dy_dt=beta*xtau*vtau-mu*y*z-delta*y;
dv_dt=b*delta*y-beta*x*v-theta*v*z-gamma*v;
dz_dt=s*z*y-c*z^2;
%% Differentiate and generate code
[fstr,derivs]=dde_sym2funcs(...
[dx_dt;dy_dt;dv_dt;dz_dt],... % n x 1 array of derivative symbolic expressions
[x,xtau;y,ytau;v,vtau;z,ztau],... % n x (ntau+1) array of symbols for states
par,... % 1 x np (or np x 1) array of symbols used for parameters
'filename','sym_virotherapy'); % argument output file