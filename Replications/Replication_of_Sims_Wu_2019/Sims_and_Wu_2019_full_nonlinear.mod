//*************************************************************************
// Dynare code to replicate Sims and Wu (2019) "The Four Equation New Keynesian Model",
//   NBER Working Paper Series 26067.
// Programmer: Camilo Marchesini
//*************************************************************************

//******************************
// Full Nonlinear model.
   (still need to add expression for output gap. Output growth in the Taylor Rule for now)
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
   pinf
   lamb
   lambb
   Rs
   Rb
   Rre
   omega
   pm
   w
   pstar
   vp
   Q
   A
   thet
   L
   C
   Cb
   Y
   b
   b_FI
   b_CB
   QE
   re
   s
   x1
   x2
  % Yf % Flex price output
  % X % Output gap 
;

// Exogenous variables 
varexo 
       eps_thet
       eps_A
       eps_QE     
;

// Parameters
parameters 
           BETA
           SIGMA
           CHI
           PSI
           BETAB
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
           RHO_R
           RHO_A
           RHO_THET
           RHO_QE
           SIGM_A
           SIGM_THET
           SIGM_QE      
;

// Parameter values 

BETA=0.995;
SIGMA=1;
CHI=1;
PSI=1.36;
BETAB=0.99;
PINF_ss=1;
EPSILONP=11;
KAPPA=1-40^(-1);
PHI=0.75;
L_ss=1;
THET_ss=5;
X_FI=0.046;
PHI_R=0.8;
PHI_PI=1.5;
PHI_X=0; % Should be the coefficient on output gap, though. Will derive the flixible price output later

RHO_R=0.8;
RHO_A=0.8;
RHO_THET=0.8;
RHO_QE=0.8;

SIGM_A=0.01;
SIGM_THET=0.01;
SIGM_QE=0.01;

model;

// Model-local variables and composite parameters ************************

#Rs_ss=1/BETA;
#Rb_ss=1/BETAB;
#pm_ss=(EPSILONP-1)/EPSILONP;
#w_ss=pm_ss;
#Q_ss=1/(Rb_ss-KAPPA);
#C_ss=((PSI*L_ss^CHI)/w_ss)^(-1/SIGMA);
#Cb_ss=1-C_ss;
#b_ss=Cb_ss/Q_ss;
#b_FI_ss=(THET_ss*X_FI)/Q_ss;
#b_CB_ss=b_ss-b_FI_ss;
#QE_ss=Q_ss*b_CB_ss;

// ************************************************************************
[name='FOC labour']
PSI*L^CHI=C^(-SIGMA)*w;

[name='Stochastic discount factor parent']
lamb=BETA*(C/C(-1))^(-SIGMA);

[name='Lucas asset pricing']
1=Rs*lamb(+1)*pinf(+1)^(-1);

[name='Long term bond pricing']
Rb=(1+KAPPA*Q)/(Q(-1));

[name='Stochastic discount factor children']
lambb=BETAB*(Cb/Cb(-1))^(-SIGMA);

[name='Lucas Asset pricing long term']
1=lambb(+1)*Rb(+1)*pinf(+1)^(-1);

[name='Spread long-short']
lamb(+1)*pinf(+1)^(-1)*(Rb(+1)-Rs)=omega;

[name='Spread reserve rate-deposit rate']
lamb(+1)*pinf(+1)^(-1)*(Rre-Rs)=0;

[name='Financial conditions']
Q*b_FI=thet*X_FI;

[name='Nominal bonds']
Q*(b_FI-KAPPA*pinf^(-1)*b_FI(-1))+re=s+X_FI;

[name='real wage']
w=pm*exp(A);

[name='Optimal reset price']
EPSILONP*x1=x2*(EPSILONP-1);

[name='Real marginal costs']
x1=pm*Y+lamb(+1)*PHI*pinf(+1)^(EPSILONP)*x1(+1);

[name='Real marginal revenues']
x2=Y+lamb(+1)*PHI*pinf(+1)^(EPSILONP-1)*x2(+1);

[name='Evolution of inflation']
1=(1-PHI)*pstar^(1-EPSILONP)+PHI*pinf^(EPSILONP-1);

[name='price dispersion']
vp=(1-PHI)*pstar^(-EPSILONP)+PHI*pinf^EPSILONP*vp(-1);

[name='Aggregate production function with price dispersion']
Y*vp=exp(A)*L;

[name='Aggregate resource constraint']
Y=C+Cb;

[name='Monetary Policy']
// It should respond to output gap here, not growth. 
Rs/Rs_ss=(Rs(-1)/Rs_ss)^RHO_R*((pinf/PINF_ss)^PHI_PI*(Y/Y(-1))^PHI_X)^(1-RHO_R);

[name='Reserves']
Q*b_CB=re;

[name='QE']
QE=Q*b_CB;

[name='Consumption of the children']
Cb=Q*b;

[name='Aggregate bond holdings']
b=b_FI+b_CB;

% X=Y/Yf; % Output gap 
//************************* Shock processes *******************************
[name='technology shock']
A=RHO_A*A(-1)+eps_A;

[name='Credit conditions shock']
thet=(1-RHO_THET)*THET_ss+RHO_THET*thet(-1)+eps_thet;

[name='QE shock']
QE=(1-RHO_QE)*QE_ss+RHO_QE*QE(-1)+eps_QE;

end;



resid(1);
steady;
check;
model_diagnostics;