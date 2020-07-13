//*************************************************************************
// Dynare code to replicate Sims and Wu (2019) "The Four Equation New Keynesian Model",
//   NBER Working Paper Series 26067.
// Programmer: Camilo Marchesini
//*************************************************************************

//*************************************************************************
// Baseline Four Equation model.
//*************************************************************************


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


// Declare relevant macros.
// Run GenerateFigures.m to generate IRFs for the various shocks and compare
// 3 and 4 equation models.

@#ifndef QE_shock
    @#define QE_shock=0
@#endif

@#ifndef RF_shock
    @#define RF_shock=0
@#endif

@#ifndef THETA_shock
    @#define THETA_shock=0
@#endif

@#ifndef MP_shock
    @#define MP_shock=0
@#endif


// Endogenous variables 
var   a_theta                (long_name='Credit conditions')
      a_rf                   (long_name='Natural rate of interest')
      a_mp                   (long_name='Policy rate shock')
      qe                     (long_name='Long bond portfolio')
      x                      (long_name='Output gap')
      pinf                   (long_name='Consumer-price inflation')
      r                      (long_name='Policy rate')
;

// Exogenous variables 
varexo eps_qe                (long_name='QE shock')
       eps_rf                (long_name='Natural interest rate shock')
       eps_theta             (long_name='Credit shock')
       eps_mp                (long_name='Policy rate shock')
;

// Parameters
parameters BETA              (long_name='Discount factor')
           ZED               (long_name='Consumption share of children in the population')
           SIGMA             (long_name='Inverse Elasticity of Substitution')
           B_FIss            (long_name='Weight on leverage in IS/NKPC curves')
           B_CBss            (long_name='Weight on QE in IS/NKPC curves')
           GAMMA             (long_name='Elaticity of inflation w.r.t. marginal cost')
           ZETA              (long_name='Elaticity of output gap w.r.t. marginal cost')
           PHI_R             (long_name='Smoothing coefficient in the interest rate rule')
           PHI_PI            (long_name='Feedback coefficient on inflation in the interest rate rule')
           PHI_X             (long_name='Feedback coefficient on output gap in the interest rate rule')
           RHO_QE            (long_name='Persistence of QE shock')
           RHO_RF            (long_name='Persistence of natural rate shock')
           RHO_THETA         (long_name='Persistence of credit shock')
           RHO_MP            (long_name='Persistence of policy shock')
           SIGM_QE           (long_name='Standard deviation of QE shock')
           SIGM_RF           (long_name='Standard deviation of natural rate shock')
           SIGM_THETA        (long_name='Standard deviation of credit shock')
           SIGM_MP           (long_name='Standard deviation of policy shock')
;

// Parameter values 

BETA=0.995;
ZED=0.33;         % Switch to zero to have the standard 3-equation NK.
SIGMA=1;
B_FIss=0.70;
B_CBss=0.30;
GAMMA=0.086;
ZETA=2;
PHI_R=0.80;
PHI_PI=1.50;
PHI_X=0;

RHO_QE=0.8;
RHO_RF=0.8;
RHO_THETA=0.8;
RHO_MP=0.8;

// declare relevant shocks

  @#if QE_shock  
SIGM_QE=0.01;
SIGM_RF=0;
SIGM_THETA=0;
SIGM_MP=0;
   @#endif
  @#if RF_shock 
SIGM_QE=0;
SIGM_RF=0.01;
SIGM_THETA=0;
SIGM_MP=0;
   @#endif
   @#if THETA_shock 
SIGM_QE=0;
SIGM_RF=0;
SIGM_THETA=0.01;
SIGM_MP=0;
    @#endif
    @#if MP_shock 
SIGM_QE=0;
SIGM_RF=0;
SIGM_THETA=0;
SIGM_MP=0.01;
   @#endif

model(linear);

// Model-local variables and composite parameters

#CHI=1;
#S_F=(SIGMA*(RHO_RF-1)*(1+CHI))/(CHI*(1-ZED));


[name='IS curve']
x = x(+1)-((1-ZED)/SIGMA)*(r-pinf(+1)-a_rf)-ZED*(B_FIss*(a_theta(+1)-a_theta)+B_CBss*(qe(+1)-qe));


[name='Phillips curve']
pinf = GAMMA*ZETA*x-((ZED*GAMMA*SIGMA)/(1-ZED))*(B_FIss*a_theta+B_CBss*qe)+BETA*pinf(+1);


[name='Interest rate rule']
r = PHI_R*r(-1) +  (1-PHI_R)*( PHI_PI*pinf + PHI_X*x )+ a_mp; 
 

 [name='QE Rule']
qe = RHO_QE*qe(-1)+eps_qe;

//************************* Shock processes *******************************

[name='Shock process natural rate of interest']
a_rf=RHO_RF*a_rf(-1)+S_F*eps_rf;

[name='Shock process credit shock']
a_theta=RHO_THETA*a_theta(-1)+eps_theta;

[name='Shock process policy rate (negative)']
 a_mp=RHO_MP*a_mp(-1)-eps_mp;

end;

resid(1);

shocks;
  var eps_qe; stderr SIGM_QE*100;
  var eps_rf; stderr SIGM_RF*100;
  var eps_theta; stderr SIGM_THETA*100;
  var eps_mp; stderr SIGM_MP*100;
end;


stoch_simul(order=1, irf=20, periods=2100, nograph);
% Save.
sims_wu_4NK=oo_.irfs;
save 'sims_wu_4NK';

% Make three equation model.
set_param_value('ZED',0); % Switch to 3-equation NK model
stoch_simul(order=1, irf=20, periods=2100, nograph);
sims_wu_3NK=oo_.irfs;
save 'sims_wu_3NK';
