function varargout=sym_virotherapy(action,varargin)
%% Automatically generated with matlabFunction
% 
%#ok<*DEFNU,*INUSD,*INUSL>

switch action
  case 'ntau'
   varargout{1}=1;
   return
  case 'tp_del'
   varargout{1}=0;
   return
  case 'maxorder'
   varargout{1}=5;
   return
  case 'directional_derivative'
   varargout{1}=0;
   return
end
ind=varargin{1};
order=varargin{2};
nout=varargin{3};
f=str2func(sprintf('sym_virotherapy_%s_%d_%d',action,ind,order));
varargout=cell(nout,1);
[varargout{:}]=f(varargin{4:end});




function [out1,out2,out3,out4] = sym_virotherapy_rhs_1_0(x,y,v,z,xtau,ytau,vtau,ztau,r,K,beta,mu,delta,b,theta,gamma,s,c,tau)
%SYM_VIROTHERAPY_RHS_1_0
%    [OUT1,OUT2,OUT3,OUT4] = SYM_VIROTHERAPY_RHS_1_0(X,Y,V,Z,XTAU,YTAU,VTAU,ZTAU,R,K,BETA,MU,DELTA,B,THETA,GAMMA,S,C,TAU)

%    This function was generated by the Symbolic Math Toolbox version 9.0.
%    13-Feb-2023 22:59:56

t2 = beta.*v.*x;
t3 = -t2;
out1 = t3-r.*x.*((x+y)./K-1.0);
if nargout > 1
    out2 = -delta.*y+beta.*vtau.*xtau-mu.*y.*z;
end
if nargout > 2
    out3 = t3-gamma.*v+b.*delta.*y-theta.*v.*z;
end
if nargout > 3
    out4 = -c.*z.^2+s.*y.*z;
end


function [out1,out2,out3,out4] = sym_virotherapy_rhs_1_1(x,y,v,z,xtau,ytau,vtau,ztau,r,K,beta,mu,delta,b,theta,gamma,s,c,tau,x_d,y_d,v_d,z_d,xtau_d,ytau_d,vtau_d,ztau_d,r_d,K_d,beta_d,mu_d,delta_d,b_d,theta_d,gamma_d,s_d,c_d,tau_d)
%SYM_VIROTHERAPY_RHS_1_1
%    [OUT1,OUT2,OUT3,OUT4] = SYM_VIROTHERAPY_RHS_1_1(X,Y,V,Z,XTAU,YTAU,VTAU,ZTAU,R,K,BETA,MU,DELTA,B,THETA,GAMMA,S,C,TAU,X_D,Y_D,V_D,Z_D,XTAU_D,YTAU_D,VTAU_D,ZTAU_D,R_D,K_d,BETA_D,MU_D,DELTA_D,B_D,THETA_D,GAMMA_D,S_D,C_D,TAU_D)

%    This function was generated by the Symbolic Math Toolbox version 9.0.
%    13-Feb-2023 22:59:57

t2 = x+y;
t3 = 1.0./K;
t4 = beta_d.*v.*x;
t5 = beta.*v.*x_d;
t6 = beta.*v_d.*x;
t7 = -t4;
t8 = -t5;
t9 = -t6;
t10 = t2.*t3;
t11 = t10-1.0;
out1 = t7+t8+t9-r.*x.*(t3.*(x_d+y_d)-K_d.*t3.*t10)-r.*t11.*x_d-r_d.*t11.*x;
if nargout > 1
    out2 = -delta_d.*y-delta.*y_d+beta_d.*vtau.*xtau+beta.*vtau.*xtau_d+beta.*vtau_d.*xtau-mu_d.*y.*z-mu.*y.*z_d-mu.*y_d.*z;
end
if nargout > 2
    out3 = t7+t8+t9-gamma_d.*v-gamma.*v_d+b.*delta_d.*y+b.*delta.*y_d+b_d.*delta.*y-theta_d.*v.*z-theta.*v.*z_d-theta.*v_d.*z;
end
if nargout > 3
    out4 = -c_d.*z.^2-c.*z.*z_d.*2.0+s.*y.*z_d+s.*y_d.*z+s_d.*y.*z;
end


function [out1,out2,out3,out4] = sym_virotherapy_rhs_1_2(x,y,v,z,xtau,ytau,vtau,ztau,r,K,beta,mu,delta,b,theta,gamma,s,c,tau,x_d,y_d,v_d,z_d,xtau_d,ytau_d,vtau_d,ztau_d,x_d_d,y_d_d,v_d_d,z_d_d,xtau_d_d,ytau_d_d,vtau_d_d,ztau_d_d,r_d,K_d,beta_d,mu_d,delta_d,b_d,theta_d,gamma_d,s_d,c_d,tau_d,r_d_d,K_d_d,beta_d_d,mu_d_d,delta_d_d,b_d_d,theta_d_d,gamma_d_d,s_d_d,c_d_d,tau_d_d)
%SYM_VIROTHERAPY_RHS_1_2
%    [OUT1,OUT2,OUT3,OUT4] = SYM_VIROTHERAPY_RHS_1_2(X,Y,V,Z,XTAU,YTAU,VTAU,ZTAU,R,K,BETA,MU,DELTA,B,THETA,GAMMA,S,C,TAU,X_D,Y_D,V_D,Z_D,XTAU_D,YTAU_D,VTAU_D,ZTAU_D,X_D_D,Y_D_D,V_D_D,Z_D_D,XTAU_D_D,YTAU_D_D,VTAU_D_D,ZTAU_D_D,R_D,K_d,BETA_D,MU_D,DELTA_D,B_D,THETA_D,GAMMA_D,S_D,C_D,TAU_D,R_D_D,K_d_d,BETA_D_D,MU_D_D,DELTA_D_D,B_D_D,THETA_D_D,GAMMA_D_D,S_D_D,C_D_D,TAU_D_D)

%    This function was generated by the Symbolic Math Toolbox version 9.0.
%    13-Feb-2023 22:59:59

t2 = x+y;
t3 = x_d+y_d;
t4 = x_d_d+y_d_d;
t5 = 1.0./K;
t7 = beta_d.*v.*x_d_d;
t8 = beta_d.*v_d_d.*x;
t9 = beta_d_d.*v.*x_d;
t10 = beta_d_d.*v_d.*x;
t11 = beta.*v_d.*x_d_d;
t12 = beta.*v_d_d.*x_d;
t6 = t5.^2;
t13 = -t7;
t14 = -t8;
t15 = -t9;
t16 = -t10;
t17 = -t11;
t18 = -t12;
t19 = t2.*t5;
t20 = t3.*t5;
t21 = t4.*t5;
t22 = K_d.*t2.*t6;
t23 = K_d_d.*t2.*t6;
t26 = t19-1.0;
t24 = -t22;
t25 = -t23;
t27 = t20+t24;
t28 = t21+t25;
out1 = t13+t14+t15+t16+t17+t18+r.*x.*(K_d.*t4.*t6+K_d_d.*t3.*t6-K_d.*K_d_d.*t5.^2.*t19.*2.0)-r.*t28.*x_d-r_d.*t28.*x-r.*t27.*x_d_d-r_d_d.*t27.*x-r_d.*t26.*x_d_d-r_d_d.*t26.*x_d;
if nargout > 1
    out2 = -delta_d.*y_d_d-delta_d_d.*y_d+beta_d.*vtau.*xtau_d_d+beta_d.*vtau_d_d.*xtau+beta_d_d.*vtau.*xtau_d+beta_d_d.*vtau_d.*xtau+beta.*vtau_d.*xtau_d_d+beta.*vtau_d_d.*xtau_d-mu_d.*y.*z_d_d-mu_d.*y_d_d.*z-mu_d_d.*y.*z_d-mu_d_d.*y_d.*z-mu.*y_d.*z_d_d-mu.*y_d_d.*z_d;
end
if nargout > 2
    out3 = t13+t14+t15+t16+t17+t18-gamma_d.*v_d_d-gamma_d_d.*v_d+b.*delta_d.*y_d_d+b.*delta_d_d.*y_d+b_d.*delta_d_d.*y+b_d_d.*delta_d.*y+b_d.*delta.*y_d_d+b_d_d.*delta.*y_d-theta_d.*v.*z_d_d-theta_d.*v_d_d.*z-theta_d_d.*v.*z_d-theta_d_d.*v_d.*z-theta.*v_d.*z_d_d-theta.*v_d_d.*z_d;
end
if nargout > 3
    out4 = c.*z_d.*z_d_d.*-2.0-c_d.*z.*z_d_d.*2.0-c_d_d.*z.*z_d.*2.0+s.*y_d.*z_d_d+s.*y_d_d.*z_d+s_d.*y.*z_d_d+s_d.*y_d_d.*z+s_d_d.*y.*z_d+s_d_d.*y_d.*z;
end


function [out1,out2,out3,out4] = sym_virotherapy_rhs_1_3(x,y,v,z,xtau,ytau,vtau,ztau,r,K,beta,mu,delta,b,theta,gamma,s,c,tau,x_d,y_d,v_d,z_d,xtau_d,ytau_d,vtau_d,ztau_d,x_d_d,y_d_d,v_d_d,z_d_d,xtau_d_d,ytau_d_d,vtau_d_d,ztau_d_d,x_d_d_d,y_d_d_d,v_d_d_d,z_d_d_d,xtau_d_d_d,ytau_d_d_d,vtau_d_d_d,ztau_d_d_d,r_d,K_d,beta_d,mu_d,delta_d,b_d,theta_d,gamma_d,s_d,c_d,tau_d,r_d_d,K_d_d,beta_d_d,mu_d_d,delta_d_d,b_d_d,theta_d_d,gamma_d_d,s_d_d,c_d_d,tau_d_d,r_d_d_d,K_d_d_d,beta_d_d_d,mu_d_d_d,delta_d_d_d,b_d_d_d,theta_d_d_d,gamma_d_d_d,s_d_d_d,c_d_d_d,tau_d_d_d)
%SYM_VIROTHERAPY_RHS_1_3
%    [OUT1,OUT2,OUT3,OUT4] = SYM_VIROTHERAPY_RHS_1_3(X,Y,V,Z,XTAU,YTAU,VTAU,ZTAU,R,K,BETA,MU,DELTA,B,THETA,GAMMA,S,C,TAU,X_D,Y_D,V_D,Z_D,XTAU_D,YTAU_D,VTAU_D,ZTAU_D,X_D_D,Y_D_D,V_D_D,Z_D_D,XTAU_D_D,YTAU_D_D,VTAU_D_D,ZTAU_D_D,X_D_D_D,Y_D_D_D,V_D_D_D,Z_D_D_D,XTAU_D_D_D,YTAU_D_D_D,VTAU_D_D_D,ZTAU_D_D_D,R_D,K_d,BETA_D,MU_D,DELTA_D,B_D,THETA_D,GAMMA_D,S_D,C_D,TAU_D,R_D_D,K_d_d,BETA_D_D,MU_D_D,DELTA_D_D,B_D_D,THETA_D_D,GAMMA_D_D,S_D_D,C_D_D,TAU_D_D,R_D_D_D,K_d_d_d,BETA_D_D_D,MU_D_D_D,DELTA_D_D_D,B_D_D_D,THETA_D_D_D,GAMMA_D_D_D,S_D_D_D,C_D_D_D,TAU_D_D_D)

%    This function was generated by the Symbolic Math Toolbox version 9.0.
%    13-Feb-2023 23:00:02

t2 = x+y;
t3 = x_d+y_d;
t4 = x_d_d+y_d_d;
t5 = x_d_d_d+y_d_d_d;
t6 = 1.0./K;
t9 = beta_d.*v_d_d.*x_d_d_d;
t10 = beta_d.*v_d_d_d.*x_d_d;
t11 = beta_d_d.*v_d.*x_d_d_d;
t12 = beta_d_d.*v_d_d_d.*x_d;
t13 = beta_d_d_d.*v_d.*x_d_d;
t14 = beta_d_d_d.*v_d_d.*x_d;
t7 = t6.^2;
t8 = t6.^3;
t15 = -t9;
t16 = -t10;
t17 = -t11;
t18 = -t12;
t19 = -t13;
t20 = -t14;
t21 = t3.*t6;
t22 = t4.*t6;
t23 = t5.*t6;
t24 = K_d.*t2.*t7;
t25 = K_d_d.*t2.*t7;
t26 = K_d_d_d.*t2.*t7;
t27 = K_d_d.*t3.*t7;
t28 = K_d.*t4.*t7;
t29 = K_d_d_d.*t3.*t7;
t30 = K_d.*t5.*t7;
t31 = K_d_d_d.*t4.*t7;
t32 = K_d_d.*t5.*t7;
t36 = K_d.*K_d_d.*t2.*t8.*2.0;
t37 = K_d.*K_d_d_d.*t2.*t8.*2.0;
t38 = K_d_d.*K_d_d_d.*t2.*t8.*2.0;
t33 = -t24;
t34 = -t25;
t35 = -t26;
t39 = -t36;
t40 = -t37;
t41 = -t38;
t42 = t21+t33;
t43 = t22+t34;
t44 = t23+t35;
t45 = t27+t28+t39;
t46 = t29+t30+t40;
t47 = t31+t32+t41;
out1 = t15+t16+t17+t18+t19+t20-r.*x.*(K_d.*K_d_d.*t5.*t8.*2.0+K_d.*K_d_d_d.*t4.*t8.*2.0+K_d_d.*K_d_d_d.*t3.*t8.*2.0-K_d_d.*K_d_d_d.*t7.*t24.*6.0)+r.*t47.*x_d+r_d.*t47.*x+r.*t46.*x_d_d+r_d_d.*t46.*x+r.*t45.*x_d_d_d+r_d_d_d.*t45.*x-r_d.*t44.*x_d_d-r_d_d.*t44.*x_d-r_d.*t43.*x_d_d_d-r_d_d_d.*t43.*x_d-r_d_d.*t42.*x_d_d_d-r_d_d_d.*t42.*x_d_d;
if nargout > 1
    out2 = beta_d.*vtau_d_d.*xtau_d_d_d+beta_d.*vtau_d_d_d.*xtau_d_d+beta_d_d.*vtau_d.*xtau_d_d_d+beta_d_d.*vtau_d_d_d.*xtau_d+beta_d_d_d.*vtau_d.*xtau_d_d+beta_d_d_d.*vtau_d_d.*xtau_d-mu_d.*y_d_d.*z_d_d_d-mu_d.*y_d_d_d.*z_d_d-mu_d_d.*y_d.*z_d_d_d-mu_d_d.*y_d_d_d.*z_d-mu_d_d_d.*y_d.*z_d_d-mu_d_d_d.*y_d_d.*z_d;
end
if nargout > 2
    out3 = t15+t16+t17+t18+t19+t20+b_d.*delta_d_d.*y_d_d_d+b_d.*delta_d_d_d.*y_d_d+b_d_d.*delta_d.*y_d_d_d+b_d_d.*delta_d_d_d.*y_d+b_d_d_d.*delta_d.*y_d_d+b_d_d_d.*delta_d_d.*y_d-theta_d.*v_d_d.*z_d_d_d-theta_d.*v_d_d_d.*z_d_d-theta_d_d.*v_d.*z_d_d_d-theta_d_d.*v_d_d_d.*z_d-theta_d_d_d.*v_d.*z_d_d-theta_d_d_d.*v_d_d.*z_d;
end
if nargout > 3
    out4 = c_d.*z_d_d.*z_d_d_d.*-2.0-c_d_d.*z_d.*z_d_d_d.*2.0-c_d_d_d.*z_d.*z_d_d.*2.0+s_d.*y_d_d.*z_d_d_d+s_d.*y_d_d_d.*z_d_d+s_d_d.*y_d.*z_d_d_d+s_d_d.*y_d_d_d.*z_d+s_d_d_d.*y_d.*z_d_d+s_d_d_d.*y_d_d.*z_d;
end


function [out1,out2,out3,out4] = sym_virotherapy_rhs_1_4(x,y,v,z,xtau,ytau,vtau,ztau,r,K,beta,mu,delta,b,theta,gamma,s,c,tau,x_d,y_d,v_d,z_d,xtau_d,ytau_d,vtau_d,ztau_d,x_d_d,y_d_d,v_d_d,z_d_d,xtau_d_d,ytau_d_d,vtau_d_d,ztau_d_d,x_d_d_d,y_d_d_d,v_d_d_d,z_d_d_d,xtau_d_d_d,ytau_d_d_d,vtau_d_d_d,ztau_d_d_d,x_d_d_d_d,y_d_d_d_d,v_d_d_d_d,z_d_d_d_d,xtau_d_d_d_d,ytau_d_d_d_d,vtau_d_d_d_d,ztau_d_d_d_d,r_d,K_d,beta_d,mu_d,delta_d,b_d,theta_d,gamma_d,s_d,c_d,tau_d,r_d_d,K_d_d,beta_d_d,mu_d_d,delta_d_d,b_d_d,theta_d_d,gamma_d_d,s_d_d,c_d_d,tau_d_d,r_d_d_d,K_d_d_d,beta_d_d_d,mu_d_d_d,delta_d_d_d,b_d_d_d,theta_d_d_d,gamma_d_d_d,s_d_d_d,c_d_d_d,tau_d_d_d,r_d_d_d_d,K_d_d_d_d,beta_d_d_d_d,mu_d_d_d_d,delta_d_d_d_d,b_d_d_d_d,theta_d_d_d_d,gamma_d_d_d_d,s_d_d_d_d,c_d_d_d_d,tau_d_d_d_d)
%SYM_VIROTHERAPY_RHS_1_4
%    [OUT1,OUT2,OUT3,OUT4] = SYM_VIROTHERAPY_RHS_1_4(X,Y,V,Z,XTAU,YTAU,VTAU,ZTAU,R,K,BETA,MU,DELTA,B,THETA,GAMMA,S,C,TAU,X_D,Y_D,V_D,Z_D,XTAU_D,YTAU_D,VTAU_D,ZTAU_D,X_D_D,Y_D_D,V_D_D,Z_D_D,XTAU_D_D,YTAU_D_D,VTAU_D_D,ZTAU_D_D,X_D_D_D,Y_D_D_D,V_D_D_D,Z_D_D_D,XTAU_D_D_D,YTAU_D_D_D,VTAU_D_D_D,ZTAU_D_D_D,X_D_D_D_D,Y_D_D_D_D,V_D_D_D_D,Z_D_D_D_D,XTAU_D_D_D_D,YTAU_D_D_D_D,VTAU_D_D_D_D,ZTAU_D_D_D_D,R_D,K_d,BETA_D,MU_D,DELTA_D,B_D,THETA_D,GAMMA_D,S_D,C_D,TAU_D,R_D_D,K_d_d,BETA_D_D,MU_D_D,DELTA_D_D,B_D_D,THETA_D_D,GAMMA_D_D,S_D_D,C_D_D,TAU_D_D,R_D_D_D,K_d_d_d,BETA_D_D_D,MU_D_D_D,DELTA_D_D_D,B_D_D_D,THETA_D_D_D,GAMMA_D_D_D,S_D_D_D,C_D_D_D,TAU_D_D_D,R_D_D_D_D,K_d_d_d_d,BETA_D_D_D_D,MU_D_D_D_D,DELTA_D_D_D_D,B_D_D_D_D,THETA_D_D_D_D,GAMMA_D_D_D_D,S_D_D_D_D,C_D_D_D_D,TAU_D_D_D_D)

%    This function was generated by the Symbolic Math Toolbox version 9.0.
%    13-Feb-2023 23:00:05

t2 = x+y;
t3 = x_d+y_d;
t4 = x_d_d+y_d_d;
t5 = x_d_d_d+y_d_d_d;
t6 = x_d_d_d_d+y_d_d_d_d;
t7 = 1.0./K.^2;
t8 = 1.0./K.^3;
t9 = t7.^2;
t10 = K_d_d.*t3.*t7;
t11 = K_d.*t4.*t7;
t12 = K_d_d_d.*t3.*t7;
t13 = K_d_d_d_d.*t3.*t7;
t14 = K_d.*t5.*t7;
t15 = K_d_d_d.*t4.*t7;
t16 = K_d_d.*t5.*t7;
t17 = K_d_d_d_d.*t4.*t7;
t18 = K_d.*t6.*t7;
t19 = K_d_d.*t6.*t7;
t20 = K_d_d_d_d.*t5.*t7;
t21 = K_d_d_d.*t6.*t7;
t22 = K_d.*K_d_d.*t2.*t8.*2.0;
t23 = K_d.*K_d_d_d.*t2.*t8.*2.0;
t24 = K_d.*K_d_d_d_d.*t2.*t8.*2.0;
t25 = K_d_d.*K_d_d_d.*t2.*t8.*2.0;
t26 = K_d_d.*K_d_d_d_d.*t2.*t8.*2.0;
t27 = K_d_d_d.*K_d_d_d_d.*t2.*t8.*2.0;
t28 = K_d_d.*K_d_d_d.*t3.*t8.*2.0;
t29 = K_d.*K_d_d_d.*t4.*t8.*2.0;
t30 = K_d_d.*K_d_d_d_d.*t3.*t8.*2.0;
t31 = K_d.*K_d_d.*t5.*t8.*2.0;
t32 = K_d.*K_d_d_d_d.*t4.*t8.*2.0;
t33 = K_d_d_d.*K_d_d_d_d.*t3.*t8.*2.0;
t34 = K_d.*K_d_d.*t6.*t8.*2.0;
t35 = K_d.*K_d_d_d_d.*t5.*t8.*2.0;
t36 = K_d_d_d.*K_d_d_d_d.*t4.*t8.*2.0;
t37 = K_d.*K_d_d_d.*t6.*t8.*2.0;
t38 = K_d_d.*K_d_d_d_d.*t5.*t8.*2.0;
t39 = K_d_d.*K_d_d_d.*t6.*t8.*2.0;
t40 = K_d.*K_d_d.*K_d_d_d.*t2.*t9.*6.0;
t41 = K_d.*K_d_d.*K_d_d_d_d.*t2.*t9.*6.0;
t42 = K_d.*K_d_d_d.*K_d_d_d_d.*t2.*t9.*6.0;
t43 = K_d_d.*K_d_d_d.*K_d_d_d_d.*t2.*t9.*6.0;
t44 = -t22;
t45 = -t23;
t46 = -t24;
t47 = -t25;
t48 = -t26;
t49 = -t27;
t50 = -t40;
t51 = -t41;
t52 = -t42;
t53 = -t43;
t54 = t10+t11+t44;
t55 = t12+t14+t45;
t56 = t13+t18+t46;
t57 = t15+t16+t47;
t58 = t17+t19+t48;
t59 = t20+t21+t49;
t60 = t28+t29+t31+t50;
t61 = t30+t32+t34+t51;
t62 = t33+t35+t37+t52;
t63 = t36+t38+t39+t53;
out1 = -r.*t63.*x_d-r_d.*t63.*x-r.*t62.*x_d_d-r_d_d.*t62.*x-r.*t61.*x_d_d_d-r_d_d_d.*t61.*x-r.*t60.*x_d_d_d_d-r_d_d_d_d.*t60.*x+r_d.*t59.*x_d_d+r_d_d.*t59.*x_d+r_d.*t58.*x_d_d_d+r_d_d_d.*t58.*x_d+r_d_d.*t56.*x_d_d_d+r_d_d_d.*t56.*x_d_d+r_d.*t57.*x_d_d_d_d+r_d_d_d_d.*t57.*x_d+r_d_d.*t55.*x_d_d_d_d+r_d_d_d_d.*t55.*x_d_d+r_d_d_d.*t54.*x_d_d_d_d+r_d_d_d_d.*t54.*x_d_d_d+r.*x.*(K_d.*K_d_d.*K_d_d_d.*t6.*t9.*6.0+K_d.*K_d_d.*K_d_d_d_d.*t5.*t9.*6.0+K_d.*K_d_d_d.*K_d_d_d_d.*t4.*t9.*6.0+K_d_d.*K_d_d_d.*K_d_d_d_d.*t3.*t9.*6.0-1.0./K.^5.*K_d.*K_d_d.*K_d_d_d.*K_d_d_d_d.*t2.*2.4e+1);
if nargout > 1
    out2 = 0.0;
end
if nargout > 2
    out3 = 0.0;
end
if nargout > 3
    out4 = 0.0;
end


function [out1,out2,out3,out4] = sym_virotherapy_rhs_1_5(x,y,v,z,xtau,ytau,vtau,ztau,r,K,beta,mu,delta,b,theta,gamma,s,c,tau,x_d,y_d,v_d,z_d,xtau_d,ytau_d,vtau_d,ztau_d,x_d_d,y_d_d,v_d_d,z_d_d,xtau_d_d,ytau_d_d,vtau_d_d,ztau_d_d,x_d_d_d,y_d_d_d,v_d_d_d,z_d_d_d,xtau_d_d_d,ytau_d_d_d,vtau_d_d_d,ztau_d_d_d,x_d_d_d_d,y_d_d_d_d,v_d_d_d_d,z_d_d_d_d,xtau_d_d_d_d,ytau_d_d_d_d,vtau_d_d_d_d,ztau_d_d_d_d,x_d_d_d_d_d,y_d_d_d_d_d,v_d_d_d_d_d,z_d_d_d_d_d,xtau_d_d_d_d_d,ytau_d_d_d_d_d,vtau_d_d_d_d_d,ztau_d_d_d_d_d,r_d,K_d,beta_d,mu_d,delta_d,b_d,theta_d,gamma_d,s_d,c_d,tau_d,r_d_d,K_d_d,beta_d_d,mu_d_d,delta_d_d,b_d_d,theta_d_d,gamma_d_d,s_d_d,c_d_d,tau_d_d,r_d_d_d,K_d_d_d,beta_d_d_d,mu_d_d_d,delta_d_d_d,b_d_d_d,theta_d_d_d,gamma_d_d_d,s_d_d_d,c_d_d_d,tau_d_d_d,r_d_d_d_d,K_d_d_d_d,beta_d_d_d_d,mu_d_d_d_d,delta_d_d_d_d,b_d_d_d_d,theta_d_d_d_d,gamma_d_d_d_d,s_d_d_d_d,c_d_d_d_d,tau_d_d_d_d,r_d_d_d_d_d,K_d_d_d_d_d,beta_d_d_d_d_d,mu_d_d_d_d_d,delta_d_d_d_d_d,b_d_d_d_d_d,theta_d_d_d_d_d,gamma_d_d_d_d_d,s_d_d_d_d_d,c_d_d_d_d_d,tau_d_d_d_d_d)
%SYM_VIROTHERAPY_RHS_1_5
%    [OUT1,OUT2,OUT3,OUT4] = SYM_VIROTHERAPY_RHS_1_5(X,Y,V,Z,XTAU,YTAU,VTAU,ZTAU,R,K,BETA,MU,DELTA,B,THETA,GAMMA,S,C,TAU,X_D,Y_D,V_D,Z_D,XTAU_D,YTAU_D,VTAU_D,ZTAU_D,X_D_D,Y_D_D,V_D_D,Z_D_D,XTAU_D_D,YTAU_D_D,VTAU_D_D,ZTAU_D_D,X_D_D_D,Y_D_D_D,V_D_D_D,Z_D_D_D,XTAU_D_D_D,YTAU_D_D_D,VTAU_D_D_D,ZTAU_D_D_D,X_D_D_D_D,Y_D_D_D_D,V_D_D_D_D,Z_D_D_D_D,XTAU_D_D_D_D,YTAU_D_D_D_D,VTAU_D_D_D_D,ZTAU_D_D_D_D,X_D_D_D_D_D,Y_D_D_D_D_D,V_D_D_D_D_D,Z_D_D_D_D_D,XTAU_D_D_D_D_D,YTAU_D_D_D_D_D,VTAU_D_D_D_D_D,ZTAU_D_D_D_D_D,R_D,K_d,BETA_D,MU_D,DELTA_D,B_D,THETA_D,GAMMA_D,S_D,C_D,TAU_D,R_D_D,K_d_d,BETA_D_D,MU_D_D,DELTA_D_D,B_D_D,THETA_D_D,GAMMA_D_D,S_D_D,C_D_D,TAU_D_D,R_D_D_D,K_d_d_d,BETA_D_D_D,MU_D_D_D,DELTA_D_D_D,B_D_D_D,THETA_D_D_D,GAMMA_D_D_D,S_D_D_D,C_D_D_D,TAU_D_D_D,R_D_D_D_D,K_d_d_d_d,BETA_D_D_D_D,MU_D_D_D_D,DELTA_D_D_D_D,B_D_D_D_D,THETA_D_D_D_D,GAMMA_D_D_D_D,S_D_D_D_D,C_D_D_D_D,TAU_D_D_D_D,R_D_D_D_D_D,K_d_d_d_d_d,BETA_D_D_D_D_D,MU_D_D_D_D_D,DELTA_D_D_D_D_D,B_D_D_D_D_D,THETA_D_D_D_D_D,GAMMA_D_D_D_D_D,S_D_D_D_D_D,C_D_D_D_D_D,TAU_D_D_D_D_D)

%    This function was generated by the Symbolic Math Toolbox version 9.0.
%    13-Feb-2023 23:00:08

t2 = x+y;
t3 = x_d+y_d;
t4 = x_d_d+y_d_d;
t5 = x_d_d_d+y_d_d_d;
t6 = x_d_d_d_d+y_d_d_d_d;
t7 = x_d_d_d_d_d+y_d_d_d_d_d;
t8 = 1.0./K.^3;
t9 = 1.0./K.^4;
t10 = 1.0./K.^5;
t11 = K_d_d.*K_d_d_d.*t3.*t8.*2.0;
t12 = K_d.*K_d_d_d.*t4.*t8.*2.0;
t13 = K_d_d.*K_d_d_d_d.*t3.*t8.*2.0;
t14 = K_d.*K_d_d.*t5.*t8.*2.0;
t15 = K_d.*K_d_d_d_d.*t4.*t8.*2.0;
t16 = K_d_d.*K_d_d_d_d_d.*t3.*t8.*2.0;
t17 = K_d_d_d.*K_d_d_d_d.*t3.*t8.*2.0;
t18 = K_d.*K_d_d_d_d_d.*t4.*t8.*2.0;
t19 = K_d_d_d.*K_d_d_d_d_d.*t3.*t8.*2.0;
t20 = K_d.*K_d_d.*t6.*t8.*2.0;
t21 = K_d.*K_d_d_d_d.*t5.*t8.*2.0;
t22 = K_d_d_d.*K_d_d_d_d.*t4.*t8.*2.0;
t23 = K_d_d_d_d.*K_d_d_d_d_d.*t3.*t8.*2.0;
t24 = K_d.*K_d_d_d.*t6.*t8.*2.0;
t25 = K_d.*K_d_d_d_d_d.*t5.*t8.*2.0;
t26 = K_d_d.*K_d_d_d_d.*t5.*t8.*2.0;
t27 = K_d_d_d.*K_d_d_d_d_d.*t4.*t8.*2.0;
t28 = K_d.*K_d_d.*t7.*t8.*2.0;
t29 = K_d_d.*K_d_d_d.*t6.*t8.*2.0;
t30 = K_d_d.*K_d_d_d_d_d.*t5.*t8.*2.0;
t31 = K_d_d_d_d.*K_d_d_d_d_d.*t4.*t8.*2.0;
t32 = K_d.*K_d_d_d.*t7.*t8.*2.0;
t33 = K_d.*K_d_d_d_d_d.*t6.*t8.*2.0;
t34 = K_d.*K_d_d_d_d.*t7.*t8.*2.0;
t35 = K_d_d.*K_d_d_d.*t7.*t8.*2.0;
t36 = K_d_d.*K_d_d_d_d_d.*t6.*t8.*2.0;
t37 = K_d_d_d_d.*K_d_d_d_d_d.*t5.*t8.*2.0;
t38 = K_d_d.*K_d_d_d_d.*t7.*t8.*2.0;
t39 = K_d_d_d.*K_d_d_d_d_d.*t6.*t8.*2.0;
t40 = K_d_d_d.*K_d_d_d_d.*t7.*t8.*2.0;
t41 = K_d.*K_d_d.*K_d_d_d.*t2.*t9.*6.0;
t42 = K_d.*K_d_d.*K_d_d_d_d.*t2.*t9.*6.0;
t43 = K_d.*K_d_d.*K_d_d_d_d_d.*t2.*t9.*6.0;
t44 = K_d.*K_d_d_d.*K_d_d_d_d.*t2.*t9.*6.0;
t45 = K_d.*K_d_d_d.*K_d_d_d_d_d.*t2.*t9.*6.0;
t46 = K_d_d.*K_d_d_d.*K_d_d_d_d.*t2.*t9.*6.0;
t47 = K_d.*K_d_d_d_d.*K_d_d_d_d_d.*t2.*t9.*6.0;
t48 = K_d_d.*K_d_d_d.*K_d_d_d_d_d.*t2.*t9.*6.0;
t49 = K_d_d.*K_d_d_d_d.*K_d_d_d_d_d.*t2.*t9.*6.0;
t50 = K_d_d_d.*K_d_d_d_d.*K_d_d_d_d_d.*t2.*t9.*6.0;
t51 = K_d_d.*K_d_d_d.*K_d_d_d_d.*t3.*t9.*6.0;
t52 = K_d.*K_d_d_d.*K_d_d_d_d.*t4.*t9.*6.0;
t53 = K_d_d.*K_d_d_d.*K_d_d_d_d_d.*t3.*t9.*6.0;
t54 = K_d.*K_d_d.*K_d_d_d_d.*t5.*t9.*6.0;
t55 = K_d.*K_d_d_d.*K_d_d_d_d_d.*t4.*t9.*6.0;
t56 = K_d_d.*K_d_d_d_d.*K_d_d_d_d_d.*t3.*t9.*6.0;
t57 = K_d.*K_d_d.*K_d_d_d.*t6.*t9.*6.0;
t58 = K_d.*K_d_d.*K_d_d_d_d_d.*t5.*t9.*6.0;
t59 = K_d.*K_d_d_d_d.*K_d_d_d_d_d.*t4.*t9.*6.0;
t60 = K_d_d_d.*K_d_d_d_d.*K_d_d_d_d_d.*t3.*t9.*6.0;
t61 = K_d.*K_d_d.*K_d_d_d.*t7.*t9.*6.0;
t62 = K_d.*K_d_d.*K_d_d_d_d_d.*t6.*t9.*6.0;
t63 = K_d.*K_d_d_d_d.*K_d_d_d_d_d.*t5.*t9.*6.0;
t64 = K_d_d_d.*K_d_d_d_d.*K_d_d_d_d_d.*t4.*t9.*6.0;
t65 = K_d.*K_d_d.*K_d_d_d_d.*t7.*t9.*6.0;
t66 = K_d.*K_d_d_d.*K_d_d_d_d_d.*t6.*t9.*6.0;
t67 = K_d_d.*K_d_d_d_d.*K_d_d_d_d_d.*t5.*t9.*6.0;
t68 = K_d.*K_d_d_d.*K_d_d_d_d.*t7.*t9.*6.0;
t69 = K_d_d.*K_d_d_d.*K_d_d_d_d_d.*t6.*t9.*6.0;
t70 = K_d_d.*K_d_d_d.*K_d_d_d_d.*t7.*t9.*6.0;
t81 = K_d.*K_d_d.*K_d_d_d.*K_d_d_d_d.*t2.*t10.*2.4e+1;
t82 = K_d.*K_d_d.*K_d_d_d.*K_d_d_d_d_d.*t2.*t10.*2.4e+1;
t83 = K_d.*K_d_d.*K_d_d_d_d.*K_d_d_d_d_d.*t2.*t10.*2.4e+1;
t84 = K_d.*K_d_d_d.*K_d_d_d_d.*K_d_d_d_d_d.*t2.*t10.*2.4e+1;
t85 = K_d_d.*K_d_d_d.*K_d_d_d_d.*K_d_d_d_d_d.*t2.*t10.*2.4e+1;
t71 = -t41;
t72 = -t42;
t73 = -t43;
t74 = -t44;
t75 = -t45;
t76 = -t46;
t77 = -t47;
t78 = -t48;
t79 = -t49;
t80 = -t50;
t86 = -t81;
t87 = -t82;
t88 = -t83;
t89 = -t84;
t90 = -t85;
t91 = t11+t12+t14+t71;
t92 = t13+t15+t20+t72;
t93 = t16+t18+t28+t73;
t94 = t17+t21+t24+t74;
t95 = t19+t25+t32+t75;
t96 = t22+t26+t29+t76;
t97 = t23+t33+t34+t77;
t98 = t27+t30+t35+t78;
t99 = t31+t36+t38+t79;
t100 = t37+t39+t40+t80;
t101 = t51+t52+t54+t57+t86;
t102 = t53+t55+t58+t61+t87;
t103 = t56+t59+t62+t65+t88;
t104 = t60+t63+t66+t68+t89;
t105 = t64+t67+t69+t70+t90;
et1 = -r.*x.*(K_d.*K_d_d.*K_d_d_d.*K_d_d_d_d.*t7.*t10.*2.4e+1+K_d.*K_d_d.*K_d_d_d.*K_d_d_d_d_d.*t6.*t10.*2.4e+1+K_d.*K_d_d.*K_d_d_d_d.*K_d_d_d_d_d.*t5.*t10.*2.4e+1+K_d.*K_d_d_d.*K_d_d_d_d.*K_d_d_d_d_d.*t4.*t10.*2.4e+1+K_d_d.*K_d_d_d.*K_d_d_d_d.*K_d_d_d_d_d.*t3.*t10.*2.4e+1-K_d.*K_d_d.*K_d_d_d.*K_d_d_d_d.*K_d_d_d_d_d.*t2.*t8.^2.*1.2e+2)+r.*t105.*x_d+r_d.*t105.*x+r.*t104.*x_d_d+r_d_d.*t104.*x+r.*t103.*x_d_d_d+r_d_d_d.*t103.*x+r.*t102.*x_d_d_d_d+r_d_d_d_d.*t102.*x+r.*t101.*x_d_d_d_d_d+r_d_d_d_d_d.*t101.*x-r_d.*t100.*x_d_d-r_d_d.*t100.*x_d-r_d.*t99.*x_d_d_d-r_d_d_d.*t99.*x_d-r_d_d.*t97.*x_d_d_d-r_d_d_d.*t97.*x_d_d-r_d.*t98.*x_d_d_d_d-r_d_d_d_d.*t98.*x_d-r_d_d.*t95.*x_d_d_d_d-r_d_d_d_d.*t95.*x_d_d-r_d.*t96.*x_d_d_d_d_d-r_d_d_d_d_d.*t96.*x_d-r_d_d_d.*t93.*x_d_d_d_d-r_d_d_d_d.*t93.*x_d_d_d-r_d_d.*t94.*x_d_d_d_d_d;
et2 = -r_d_d_d_d_d.*t94.*x_d_d-r_d_d_d.*t92.*x_d_d_d_d_d-r_d_d_d_d_d.*t92.*x_d_d_d-r_d_d_d_d.*t91.*x_d_d_d_d_d-r_d_d_d_d_d.*t91.*x_d_d_d_d;
out1 = et1+et2;
if nargout > 1
    out2 = 0.0;
end
if nargout > 2
    out3 = 0.0;
end
if nargout > 3
    out4 = 0.0;
end

