//*************************************************************************
// Dynare code to replicate Sims and Wu (2019) "The Four Equation New Keynesian Model",
//   NBER Working Paper Series 26067.
// Programmer: Camilo Marchesini
//*************************************************************************

//******************************
// Full Linearised Model.
//******************************

% This implementation was written by Camilo Marchesini.

% MATLAB_R2019a and subsequent distributions. Backward compatibility untested.

%     Copyright (C) 2020  Camilo Marchesini
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.




// Endogenous variables 
var   
       l
       c
       w
       lamb
       rs
       pinf
       lambb
       rb
       q
       b_FI
       thet
       re
       s
       omega
       rre
       pstar
       x1
       x2
       pm
       y
       a
       cb
       vp
       b_CB
       b
       qe
       y_flex
       xgap       
;


// Exogenous variables 
varexo 
       eps_thet
       eps_A
       eps_QE     
;


// Parameters
parameters BETA              % (long name="Discount factor")
           ZED               % (long name="Consumption share of children in the population")
           SIGMA             % (long name="Inverse Elasticity of Substitution")
           CHI
           BETAB
           PSI
           PINF_ss
           L_ss
           THET_ss
           EPSILONP
           KAPPA
           PHI
           X_FI
           PHI_R
           PHI_PI
           PHI_X
           RHO_A
           RHO_THET
           RHO_QE
         %  SIGM_A
         %  SIGM_THET
         %  SIGM_QE      
;

// Parameter values 

BETA=0.995;
ZED=0.33;         % Switch to zero to have the standard 3-equation NK.
SIGMA=1;
BETAB=0.99;
CHI=1;
PSI=1.36;
PINF_ss=1;
EPSILONP=11;
KAPPA=1-40^(-1);
PHI=0.75;
L_ss=1;
THET_ss=5;
X_FI=0.046;
PHI_R=0.8;
PHI_PI=1.5;
PHI_X=0; 

RHO_A=0.8;
RHO_THET=0.8;
RHO_QE=0.8;

% SIGM_A=0.01;
% SIGM_THET=0.01;
% SIGM_QE=0.01;

model(linear);
//************************** Model-local variables and composite parameters

#lamb_ss=BETA;
#lambb_ss=BETAB;
#Rs_ss=1/BETA;
#Rb_ss=1/BETAB;
#Rre_ss=Rs_ss;
#omega_ss=BETA*(Rb_ss-Rs_ss);
#pm_ss=(EPSILONP-1)/EPSILONP;
#w_ss=pm_ss;
#pstar_ss=1;
#vp_ss=1;
#Q_ss=1/(Rb_ss-KAPPA);
#C_ss=((PSI*L_ss^CHI)/w_ss)^(-1/SIGMA);
#Cb_ss=1-C_ss;
#Y_ss=C_ss+Cb_ss;
#b_ss=Cb_ss/Q_ss;
#b_FI_ss=(THET_ss*X_FI)/Q_ss;
#b_CB_ss=b_ss-b_FI_ss;
#QE_ss=Q_ss*b_CB_ss;
#re_ss=QE_ss;
#s_ss=X_FI*(THET_ss*(1-KAPPA)-1)+re_ss;
#x1_ss=(Y_ss*pm_ss)/(1-BETA*PHI);
#x2_ss=Y_ss/(1-BETA*PHI);

//*************************************************************************

// pp. 34, in order.
[name='FOC labour']
CHI*l=-SIGMA*c+w;

[name='Stochastic discount factor parent']
lamb=-SIGMA*(c-c(-1));

[name='Lucas asset pricing']
lamb(+1)+rs-pinf(+1)=0;

[name='Long term bond pricing']
lambb=-SIGMA*(cb-cb(-1));

[name='Long term bond pricing']
rb=(KAPPA/Rb_ss)*q-q(-1);

[name='Stochastic discount factor children']
lambb(+1)+rb(+1)-pinf(+1)=0;

[name='Financial conditions']
q+b_FI=thet;

[name='Spread long-short']
Q_ss*b_FI_ss*(1-KAPPA*PINF_ss^(-1))*q+Q_ss*b_FI_ss*b_FI-KAPPA*PINF_ss^(-1)*Q_ss*b_FI_ss*b_FI(-1)+KAPPA*PINF_ss^(-2)*Q_ss*b_FI_ss*pinf+re*re_ss=s_ss*s;

[name='Spread long-short']
lamb(+1)-pinf(+1)+(Rb_ss/s_ss)*rb(+1)-(Rs_ss/s_ss)*rs=omega;

[name='Spread reserve rate-deposit rate']
rre=rs;

[name='Optimal reset price']
pstar=x1-x2;

[name='Real marginal costs']
x1=(1-PHI*BETA)*pm+(1-PHI*BETA)*y+PHI*BETA*lamb(+1)+EPSILONP*PHI*BETA*pinf(+1)+PHI*BETA*x1(+1);

[name='real marginal revenues']
x2=(1-PHI*BETA)*y+PHI*BETA*lamb(+1)+(EPSILONP-1)*PHI*BETA*pinf(+1)+PHI*BETA*x2(+1);

[name='real wage']
w=pm+a;

[name='Aggregate resource constraint']
y=(1-ZED)*c+ZED*cb;

[name='Aggregate production function with price dispersion']
vp+y=a+l;

[name='price dispersion']
vp=0;

[name='Evolution of inflation']
pinf=((1-PHI)/PHI)*pstar;

[name='Reserves']
q+b_CB=re;

[name='Aggregate bond holdings']
b=(b_FI_ss/b_ss)*b_FI+(b_CB_ss/b_ss)*b_CB;

[name='Consumption of the children']
cb=q+b;

[name='Monetary Policy']
rre=PHI_R*rre(-1)+(1-PHI_R)*(PHI_PI*pinf+PHI_X*xgap);

[name='QE']
qe=re;

[name='Potential Output']
y_flex=(((1+CHI)*(1-ZED))/(CHI*(1-ZED)+SIGMA))*a;

[name='Output Gap']
xgap=y-y_flex;


//************************* Shock processes *******************************
[name='QE shock']
qe=RHO_QE*qe(-1)+eps_QE;

[name='Technology Shock']
a=RHO_A*a(-1)+eps_A;

[name='Financial conditions shock']
thet=RHO_THET*thet(-1)+eps_thet;

end;

resid(1);
