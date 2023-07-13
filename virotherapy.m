%% Virutherapy model
%
% Demo illustrating how to branch off a Hopf point, branch off period
% doubling bifurcations and continue Hopf bifurcation curve in two
% parameters. Furthermore, the stability of the various branches
% is computed.
%
%% Differential equations
%
% From: Chapter 4 Bruno thesis
%
% $$\frac{dx}{dt}=rx\left(1-\frac{x+y}{K}\right)-\beta xv,$$
% $$\frac{dy}{dt}=\beta x(t-\tau)v(t-\tau)-\mu yz-\delta y,$$
% $$\frac{dv}{dt}=b\delta y-\beta xv-\theta vz-\gamma v,$$
% $$\frac{dz}{dt}=szy-cz^2.$$
%
% Parameters are (in this order) |r|, |K|, |beta|, |mu|, |delta|, |b|, |theta|, |gamma|, |s|, |c|, |tau|
%% load DDE-BifTool into MATLAB path
clear
close all
ddebiftoolpath='../../';
addpath(strcat(ddebiftoolpath,'ddebiftool'),...
strcat(ddebiftoolpath,'ddebiftool_extra_psol'),...
strcat(ddebiftoolpath,'ddebiftool_extra_nmfm'),...
strcat(ddebiftoolpath,'ddebiftool_utilities'));
format compact
set(groot, 'defaultTextInterpreter', 'LaTeX');
%
%% Initial parameters and state
parnames={'r','K','beta','mu','delta','b','theta','gamma','s','c','tau'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
bounds={'max_bound',[ind.b 550;ind.s 0.003;ind.tau 20],'max_step',[0,0.3],...
'min_bound',[ind.b 450;ind.s 0.001;ind.tau 0]};
%% Set user-defined functions
% use the right-hand side and derivatives created via symbolic toolbox
funcs=set_symfuncs(@sym_virotherapy,'sys_tau',@()ind.tau);
%
%% Construct steady-state point
r = 0.2062; K = 3000; beta = 0.00002; mu = 1/48; delta = 1/18; b = 500; theta = 2e-6; 
gamma = 0.025; s = 0.0027; c = 0.001;
stst=dde_stst_create('x',[26.8668; 9.3825; 10185.4232; 25.3327]);
stst.parameter(ind.r) = r;
stst.parameter(ind.K) = K;
stst.parameter(ind.beta) = beta;
stst.parameter(ind.mu) = mu;
stst.parameter(ind.delta) = delta;
stst.parameter(ind.b) = b;
stst.parameter(ind.theta) = theta;
stst.parameter(ind.gamma) = gamma;
stst.parameter(ind.s) = s;
stst.parameter(ind.c) = c;
stst.parameter(ind.tau) = 2;
% Compute stability
method_stst=df_mthod(funcs,'stst');
method_stst.stability.minimal_real_part=-1;
stst.stability=p_stabil(funcs,stst,method_stst.stability);
% Plot eigenvalues
figure(1); clf;
p_splot(stst); % plot its stability:
%plot(stst.stability.l1,'*')
title('Stability plot of stst')
xlabel('$\Re(\lambda)$')
ylabel('$\Im(\lambda)$')
grid on
%
stst.stability.l1
%% Initialization of branch of steady-states for b
contpar=ind.b;
steadystate_brb=SetupStst(funcs,'x',stst.x,'parameter',stst.parameter,...
'step',0.1,'contpar',contpar,'max_step',[contpar,0.3],bounds{:});
%
%% Continue steady-state branch for b
figure(2);clf;ax2=gca;
n_steps=200;
steadystate_brb=br_contn(funcs,steadystate_brb,n_steps,'plotaxis',ax2);
steadystate_brb=br_rvers(steadystate_brb);
steadystate_brb=br_contn(funcs,steadystate_brb,n_steps,'plotaxis',ax2);
xlabel('$b$')
ylabel('$x$')
grid on
%
%% Initialization of branch of steady-states for s
contpar=ind.s;
steadystate_brs=SetupStst(funcs,'x',stst.x,'parameter',stst.parameter,...
'step',0.1,'contpar',contpar,'max_step',[contpar,0.3],bounds{:});
%
%% Continue steady-state branch for s
figure(3);clf;ax2=gca;
n_steps=200;
steadystate_brs=br_contn(funcs,steadystate_brs,n_steps,'plotaxis',ax2);
steadystate_brs=br_rvers(steadystate_brs);
steadystate_brs=br_contn(funcs,steadystate_brs,n_steps,'plotaxis',ax2);
xlabel('$s$')
ylabel('$x$')
grid on
%
%% Initialization of branch of steady-states for tau
contpar=ind.tau;
steadystate_br=SetupStst(funcs,'x',stst.x,'parameter',stst.parameter,...
'step',0.1,'contpar',contpar,'max_step',[contpar,0.3],bounds{:});
%
%% Continue steady-state branch for tau
figure(4);clf;ax2=gca;
n_steps=200;
steadystate_br=br_contn(funcs,steadystate_br,n_steps,'plotaxis',ax2);
steadystate_br=br_rvers(steadystate_br);
steadystate_br=br_contn(funcs,steadystate_br,n_steps,'plotaxis',ax2);
xlabel('$\tau$')
ylabel('$x$')
grid on
%
%% Detect bifurcations on steady-state branch
[steadystate_br,~,ind_hopf,bif1types]=LocateSpecialPoints(funcs,...
steadystate_br);
nunst_eqs=GetStability(steadystate_br);
fprintf('Hopf bifurcation near point %d\n',ind_hopf);
%
steadystate_br.point(ind_hopf)
steadystate_br.point(ind_hopf).x
steadystate_br.point(ind_hopf).parameter
steadystate_br.point(ind_hopf).parameter(11)
