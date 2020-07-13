//*************************************************************************
// Dynare code to replicate Sims and Wu (2019) "The Four Equation New Keynesian Model",
//   NBER Working Paper Series 26067.
// Programmer: Camilo Marchesini
//*************************************************************************

// Baseline Four Equation model with implementation of forward guidance shocks.


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
var   a_theta                 (long name='Credit conditions')
      a_rf                    (long name='Natural rate of interest')
      a_mp                    (long name='Policy rate shock')

      a_fg_1h                 (long name='Anticipated monetary policy shock one period ahead')
      a_fg_2h                 (long name='Anticipated monetary policy shock two periods ahead')
      a_fg_3h                 (long name='Anticipated monetary policy shock three periods ahead')
      a_fg_4h                 (long name='Anticipated monetary policy shock four periods ahead')

      qe                      (long name="Long bond portfolio")
      x                       (long name="output gap")
      pinf                    (long name="Consumer-price inflation")
      r                       (long name="Policy rate")
;

// Exogenous variables 
varexo eps_qe                 (long name='QE shock')
       eps_rf                 (long name='natural interest rate shock')
       eps_theta              (long name='credit shock')
       eps_mp                 (long name='policy rate shock')

       eps_fg_1h              (long_name='Innovation FG shock one period ahead')
       eps_fg_2h              (long_name='Innovation FG shock two periods ahead')
       eps_fg_3h              (long_name='Innovation FG shock three periods ahead')
       eps_fg_4h              (long_name='Innovation FG shock four periods ahead')

;

// Parameters
parameters BETA               (long name='Discount factor')
           ZED                (long name='Consumption share of children in the population')
           SIGMA              (long name='Inverse elasticity of substitution')
           B_FIss             (long name='Weight on leverage in IS/NKPC curves')
           B_CBss             (long name='Weight on QE in IS/NKPC curves')
           GAMMA              (long name='Elaticity of inflation w.r.t. marginal cost')
           ZETA               (long name='Elaticity of output gap w.r.t. marginal cost')
           PHI_R              (long name='Smoothing coefficient in the interest rate rule')
           PHI_PI             (long name='Feedback coefficient on inflation in the interest rate rule')
           PHI_X              (long name='Feedback coefficient on output gap in the interest rate rule')
           RHO_QE             (long name='Persistence of QE shock')
           RHO_RF             (long name='Persistence of natural rate shock')
           RHO_THETA          (long name='Persistence of credit shock')
           RHO_MP             (long name='Persistence of policy shock')

           RHO_FG_1H          (long name='Persistence of FG shock one period ahead')
           RHO_FG_2H          (long name='Persistence of FG shock two periods ahead')
           RHO_FG_3H          (long name='Persistence of FG shock three periods ahead')
           RHO_FG_4H          (long name='Persistence of FG shock three periods ahead')
           
           SIGM_QE            (long name='Standard deviation of QE shock')
           SIGM_RF            (long name='Standard deviation of natural rate shock')
           SIGM_THETA         (long name='Standard deviation of credit shock')
           SIGM_MP            (long name='Standard deviation of policy shock')

           SIGM_FG_1H         (long name='Standard deviation of FG shock one period ahead')
           SIGM_FG_2H         (long name='Standard deviation of FG shock two periods ahead')
           SIGM_FG_3H         (long name='Persistence of FG shock three periods ahead')
           SIGM_FG_4H         (long name='Standard deviation of FG shock four periods ahead')
 
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

RHO_FG_1H=0.2;
RHO_FG_2H=0.2;
RHO_FG_3H=0.2;
RHO_FG_4H=0.2;

SIGM_QE=0.1;
SIGM_RF=0.1;
SIGM_THETA=0.1;
SIGM_MP=0.1;

SIGM_FG_1H=0.1;
SIGM_FG_2H=0.1;
SIGM_FG_3H=0.1;
SIGM_FG_4H=0.1;


model(linear);

// Model-local variables and composite parameters

#CHI=1;
#S_F=(SIGMA*(RHO_RF-1)*(1+CHI))/(CHI*(1-ZED));

// Eq. 2.1 
 [name='IS curve']
x = x(+1)-((1-ZED)/SIGMA)*(r-pinf(+1)-a_rf)-ZED*(B_FIss*(a_theta(+1)-a_theta)+B_CBss*(qe(+1)-qe));

// Eq. 2.2
 [name='Phillips curve']
pinf = GAMMA*ZETA*x-((ZED*GAMMA*SIGMA)/(1-ZED))*(B_FIss*a_theta+B_CBss*qe)+BETA*pinf(+1);

//Eq. 2.33
  [name='Interest rate rule']
r = PHI_R*r(-1) +  (1-PHI_R)*( PHI_PI*pinf + PHI_X*x )

+ a_mp

+ a_fg_1h(-1)
+ a_fg_2h(-2)
+ a_fg_3h(-3)
+ a_fg_4h(-4); 
 
//Eq. 2.34
 % [name='QE Rule']
qe = RHO_QE*qe(-1)+SIGM_QE*eps_qe;

//************************* Shock processes *******************************

[name='Shock process natural rate of interest']
a_rf=RHO_RF*a_rf(-1)+SIGM_RF*eps_rf;

[name='Shock process credit shock']
a_theta=RHO_THETA*a_theta(-1)+SIGM_THETA*eps_theta;

[name='Shock process policy rate (negative)']
a_mp=RHO_MP*a_mp(-1)-SIGM_MP*eps_mp;

[name='Shock process FG shock one period ahead']
a_fg_1h       =   RHO_FG_1H*a_fg_1h(-1) + SIGM_FG_1H*eps_fg_1h;

[name='Shock process FG shock two periods ahead']
a_fg_2h       =   RHO_FG_2H*a_fg_2h(-2) + SIGM_FG_2H*eps_fg_2h;

[name='Shock process FG shock three periods ahead']
a_fg_3h       =   RHO_FG_3H*a_fg_3h(-3) + SIGM_FG_3H*eps_fg_3h;

[name='Shock process FG shock four periods ahead']
a_fg_4h       =   RHO_FG_4H*a_fg_4h(-4) + SIGM_FG_4H*eps_fg_4h;

end;



shocks;
  var eps_qe;    stderr 1;
  var eps_rf;    stderr 1;
  var eps_theta; stderr 1;
  var eps_mp;    stderr 1;  
  var eps_fg_1h; stderr 1;
  var eps_fg_2h; stderr 1;
  var eps_fg_3h; stderr 1;
  var eps_fg_4h; stderr 1;
end;


resid(1);
steady;
check;
model_diagnostics;

stoch_simul(order=1, irf=20, periods=2100);