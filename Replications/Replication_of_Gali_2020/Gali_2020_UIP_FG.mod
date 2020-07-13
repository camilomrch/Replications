//*************************************************************************
// Gali, J.(2020) 
// "Uncovered Interest Parity, Forward Guidance, and the Exchange Rate",
// NBER Working Paper Series 26797.
// Programmer: Camilo Marchesini
//*************************************************************************


var 

pinf               (long_name='CPI inflation')
pinf_h             (long_name='Domestic inflation')
ygap               (long_name='Output gap')
q                  (long_name='Real exchange rate')
c                  (long_name='Consumption')
nir                (long_name='Nominal interest rate')


%e                  (long_name='Nominal exchange rate')
%plevel             (long_name='Consumer price level')
%plevel_h           (long_name='Domestic price level')
; 


varexo 
        eps_fg_1    (long_name='Anticipated monetary policy shock one quarter ahead')
        eps_fg_2    (long_name='Anticipated monetary policy shock two quarters ahead')
        eps_fg_4    (long_name='Anticipated monetary policy shock four quarters ahead')
 ;



parameters

BETA                 (long_name='Discount factor')
THETA                (long_name='Calvo probability')
SIGMA                (long_name='Coefficient of relative risk aversion')
PHI                  (long_name='Coefficient of relative risk aversion')
UPSILON              (long_name='Index of openness')
ETA                  (long_name='Elasticity of substitution between domestic and foreign goods')
PHI_PI               (long_name='Feedback coefficient on domestic inflation in the monetary policy rule')
;


% Parameter values

BETA=0.99;
THETA=0.75;
SIGMA=1;
PHI=1;
UPSILON=0.4;
ETA=2;
PHI_PI=1.75;



model(linear); 

% Composite parameters.
#LAMBDA=((1-THETA)*(1-BETA*THETA))/THETA;
#KAPPA=LAMBDA*(SIGMA+PHI);
#OMEGA=(LAMBDA*(SIGMA*ETA-1)*UPSILON*(2-UPSILON))/(1-UPSILON);
#VARTHETA=ETA*UPSILON*(1+(1/(1-UPSILON)));


[name='New Keynesian Phillips curve in a SOE']
pinf_h=BETA*pinf(+1)+KAPPA*ygap-OMEGA*q;

[name='Market clearing']
ygap=(1-UPSILON)*c+VARTHETA*q;

[name='Euler equation']
c=c(+1)-(1/SIGMA)*(nir-pinf(+1));

[name='International risk-sharing condition']
c=(1/SIGMA)*q;


[name='Relationship between CPI inflation and domestic inflation']
 pinf=pinf_h+(UPSILON/(1-UPSILON))*(q-q(-1));

[name='Constant rate until shock']
nir =  eps_fg_1 + eps_fg_2 + eps_fg_4;
% nir = PHI_PI*pinf_h + eps_fg_1 + eps_fg_2 + eps_fg_4;


%[name='Real uncovered interest parity condition']
%q=q(+1)-(nir-pinf(+1));

%[name='Domestic price level']
%plevel_h   = plevel_h(-1) + pinf_h;

%[name='CPI inflation']
%plevel     =  plevel(-1) + pinf;

%[name='Nominal exchange rate']
%e=q+plevel;

end;

model_diagnostics;

% Define anticipated monetary policy shock. 

@#ifdef FG_1 
shocks;
var eps_fg_1;
periods 1 2;
values  0 0.25; 
end;

perfect_foresight_setup(periods=20);
perfect_foresight_solver;
one_ahead=oo_.endo_simul;
@#endif

@#ifdef FG_2
shocks;
var eps_fg_2;
periods 1 2 3;
values  0 0 0.25;
end;

perfect_foresight_setup(periods=20);
perfect_foresight_solver;
two_ahead=oo_.endo_simul;
@#endif

@#ifdef FG_4
shocks; 
var eps_fg_4;
periods 1 2 3 4 5;
values  0 0 0 0 0.25; 
end;

perfect_foresight_setup(periods=20);
perfect_foresight_solver;
four_ahead=oo_.endo_simul;
@#endif
