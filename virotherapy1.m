%% Construct Hopf point manually
stst.parameter(ind.tau) = 3.7286;
stst.x = [26.8668; 9.3825; 10185.4232; 25.3327];
% Recompute stability
method_stst=df_mthod(funcs,'stst');
stst.stability=p_stabil(funcs,stst,method_stst.stability);
% Convert steady-state to Hopf point and compute the normal form
hopf=p_tohopf(funcs,stst);
hopf=nmfm_hopf(funcs,hopf);
hopf.nmfm
%
%% Branch off at Hopf bifurcation
disp('Branch off at Hopf bifurcation');
fprintf('Initial correction of periodic orbits at Hopf:\n');
[per1,suc]=SetupPsol(funcs,steadystate_br,ind_hopf,...
'print_residual_info',1,'intervals',20,'degree',4,...
'max_bound',[contpar,20],'max_step',[contpar,0.5],'matrix','full');
figure(5);clf;ax3=gca;
xlabel('$\tau$')
ylabel('Amplitude')
per1=br_contn(funcs,per1,80,'plotaxis',ax3);
grid on
%
%% Manually construct periodic solution near Hopf point
radius=0.1;
degree=4;
int_nr=20;
[psol,stpcond]=p_topsol(funcs,hopf,radius,degree,int_nr);
method=df_mthod('psol');
[psol,succ]=p_correc(funcs,psol,ind.tau,stpcond,method.point)
%% Compute stability of periodic solution
psol.stability=p_stabil(funcs,psol,method.stability);
%
psol.stability.mu(1:6)
figure(6);clf
p_splot(psol)
xlabel('$\Re(\mu)$')
ylabel('$\Im(\mu)$')
title('Floquet multipliers for periodic solution')
grid on
%
%% Manual construction and continuation of psol branch
per_orb=df_brnch(funcs,ind.tau,'psol'); % empty branch:
per_orb.parameter.max_bound=[ind.tau 20];
per_orb.parameter.max_step=[ind.tau 0.3];
deg_psol=p_topsol(funcs,hopf,0,degree,int_nr);
per_orb.point=deg_psol;
per_orb.point(2)=psol;
% compute periodic solutions branch
figure(7);clf;ax3=gca;
per_orb=br_contn(funcs,per_orb,80,'plotaxis',ax3);
xlabel('$\tau$');
ylabel('amplitude');
grid on;
%
%% Continue Hopf bifurcation in two parameters tau and b
[hbranch,suc]=SetupHopf(funcs,steadystate_br,ind_hopf,...
'contpar',[ind.tau,ind.b],'dir',ind.tau,'step',1e-1,bounds{:});
figure(8);clf;ax2=gca;
hbranch=br_contn(funcs,hbranch,200,'plotaxis',ax2);
hbranch=br_rvers(hbranch);
hbranch=br_contn(funcs,hbranch,200,'plotaxis',ax2);
xlabel('$\tau$');
ylabel('b');
grid on;
%% Compute L1 coefficient
% to find if Hopf bifurcation is supercritical (L1<0) or subcritical (L1>0)
[hbranch,hopftests,hc2_indices,hc2_types]=...
LocateSpecialPoints(funcs,hbranch);
%% Continue Hopf bifurcation in two parameters tau and s
[hbranch1,suc]=SetupHopf(funcs,steadystate_br,ind_hopf,...
'contpar',[ind.tau,ind.s],'dir',ind.tau,'step',1e-1,bounds{:});
figure(9);clf;ax2=gca;
hbranch1=br_contn(funcs,hbranch1,200,'plotaxis',ax2);
hbranch1=br_rvers(hbranch1);
hbranch1=br_contn(funcs,hbranch1,200,'plotaxis',ax2);
xlabel('$\tau$');
ylabel('s');
grid on;
%% Compute L1 coefficient
% to find if Hopf bifurcation is supercritical (L1<0) or subcritical (L1>0)
[hbranch1,hopftests,hc2_indices,hc2_types]=...
LocateSpecialPoints(funcs,hbranch1);
