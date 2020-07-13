% This script plots Figure 2 in Gali 2020 "Uncovered Interest Parity, Forward Guidance, and the Exchange Rate"

% To do: automate.


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



% Define model name.
model_name='Gali_2020_UIP_FG';

% Define macros.
macros_fg1=' -DFG_1';
macros_fg2=' -DFG_2';
macros_fg4=' -DFG_4';

% Invoke Dynare.
eval(['dynare ',model_name,macros_fg1]);
save('one_ahead','one_ahead')
% clear -model_name -macros_fg1 -macros_fg2 -macros_fg4
eval(['dynare ',model_name,macros_fg2]);
save('two_ahead','two_ahead')
% clear -model_name -macros_fg1 -macros_fg2 -macros_fg4
eval(['dynare ',model_name,macros_fg4]);
save('four_ahead','four_ahead')
% clear -model_name -macros_fg1 -macros_fg2 -macros_fg4
load('one_ahead');
load('two_ahead');
load('four_ahead');



% Length.
L=0:12;
% To add to maximum lead.
add=13;
% Number of IRF periods.
IRF_periods=L; nIRF_periods=length(IRF_periods);
% Legend names.
legend_names={'T=1','T=2','T=4'}; 

figure

subplot(2,3,1)
plot(IRF_periods,one_ahead(strmatch('nir',M_.endo_names,'exact'),M_.maximum_lag+1:M_.maximum_lag+add),'r-o'); hold on;
plot(IRF_periods,two_ahead(strmatch('nir',M_.endo_names,'exact'),M_.maximum_lag+1:M_.maximum_lag+add),'b-o'); hold on;
plot(IRF_periods,four_ahead(strmatch('nir',M_.endo_names,'exact'),M_.maximum_lag+1:M_.maximum_lag+add),'k-o'); hold on;
plot(zeros(nIRF_periods,1),'k--','HandleVisibility','off','LineWidth',1); hold off;
title('Nominal Interest rate')
xlim([0 12])

subplot(2,3,2)
plot(IRF_periods,one_ahead(strmatch('pinf',M_.endo_names,'exact'),M_.maximum_lag+1:M_.maximum_lag+add),'r-o'); hold on;
plot(IRF_periods,two_ahead(strmatch('pinf',M_.endo_names,'exact'),M_.maximum_lag+1:M_.maximum_lag+add),'b-o'); hold on;
plot(IRF_periods,four_ahead(strmatch('pinf',M_.endo_names,'exact'),M_.maximum_lag+1:M_.maximum_lag+add),'k-o'); hold on;
plot(zeros(nIRF_periods,1),'k--','HandleVisibility','off','LineWidth',1); hold off;
title('CPI Inflation')
xlim([0 12])

subplot(2,3,3)
plot(IRF_periods,one_ahead(strmatch('pinf_h',M_.endo_names,'exact'),M_.maximum_lag+1:M_.maximum_lag+add),'r-o'); hold on;
plot(IRF_periods,two_ahead(strmatch('pinf_h',M_.endo_names,'exact'),M_.maximum_lag+1:M_.maximum_lag+add),'b-o'); hold on;
plot(IRF_periods,four_ahead(strmatch('pinf_h',M_.endo_names,'exact'),M_.maximum_lag+1:M_.maximum_lag+add),'k-o'); hold on;
plot(zeros(nIRF_periods,1),'k--','HandleVisibility','off','LineWidth',1); hold off;
title('Domestic inflation')
xlim([0 12])

subplot(2,3,4)
plot(IRF_periods,one_ahead(strmatch('ygap',M_.endo_names,'exact'),M_.maximum_lag+1:M_.maximum_lag+add),'r-o'); hold on;
plot(IRF_periods,two_ahead(strmatch('ygap',M_.endo_names,'exact'),M_.maximum_lag+1:M_.maximum_lag+add),'b-o'); hold on;
plot(IRF_periods,four_ahead(strmatch('ygap',M_.endo_names,'exact'),M_.maximum_lag+1:M_.maximum_lag+add),'k-o'); hold on;
plot(zeros(nIRF_periods,1),'k--','HandleVisibility','off','LineWidth',1); hold off;
title('Output gap')
xlim([0 12])

subplot(2,3,5)
plot(IRF_periods,one_ahead(strmatch('q',M_.endo_names,'exact'),M_.maximum_lag+1:M_.maximum_lag+add),'r-o'); hold on;
plot(IRF_periods,two_ahead(strmatch('q',M_.endo_names,'exact'),M_.maximum_lag+1:M_.maximum_lag+add),'b-o'); hold on;
plot(IRF_periods,four_ahead(strmatch('q',M_.endo_names,'exact'),M_.maximum_lag+1:M_.maximum_lag+add),'k-o'); hold on;
plot(zeros(nIRF_periods,1),'k--','HandleVisibility','off','LineWidth',1); hold off;
title('Real exchange rate')
xlim([0 12])

subplot(2,3,6)
plot(IRF_periods,one_ahead(strmatch('c',M_.endo_names,'exact'),M_.maximum_lag+1:M_.maximum_lag+add),'r-o'); hold on;
plot(IRF_periods,two_ahead(strmatch('c',M_.endo_names,'exact'),M_.maximum_lag+1:M_.maximum_lag+add),'b-o'); hold on;
plot(IRF_periods,four_ahead(strmatch('c',M_.endo_names,'exact'),M_.maximum_lag+1:M_.maximum_lag+add),'k-o'); hold on;
plot(zeros(nIRF_periods,1),'k--','HandleVisibility','off','LineWidth',1); hold off;
title('Consumption')
xlim([0 12])

% Legend.
hleg = legend(legend_names{:},'Orientation','vertical');
set(hleg,...
'Location','best',...
'EdgeColor',[1 1 1],'Color','None','Box','on','FontSize',12);

% Clean-up
for d = dir(pwd).'
  if(~d.isdir && ~any(strcmp(d.name,{'Gali_2020_UIP_FG.mod','run.m','README.md'})))
    delete(fullfile(pwd, d.name))
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DO NOT EDIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p,n]=numSubplots(n)
% function [p,n]=numSubplots(n)
%
% Purpose
% Calculate how many rows and columns of sub-plots are needed to
% neatly display n subplots. 
%
% Inputs
% n - the desired number of subplots.     
%  
% Outputs
% p - a vector length 2 defining the number of rows and number of
%     columns required to show n plots.     
% [ n - the current number of subplots. This output is used only by
%       this function for a recursive call.]
%
%
%
% Example: neatly lay out 13 sub-plots
% >> p=numSubplots(13)
% p = 
%     3   5
% for i=1:13; subplot(p(1),p(2),i), pcolor(rand(10)), end 
%
%
% Rob Campbell - January 2010
   
    
while isprime(n) && n>4
    n=n+1;
end
p=factor(n);
if length(p)==1
    p=[1,p];
    return
end
while length(p)>2
    if length(p)>=4
        p(1)=p(1)*p(end-1);
        p(2)=p(2)*p(end);
        p(end-1:end)=[];
    else
        p(1)=p(1)*p(2);
        p(2)=[];
    end    
    p=sort(p);
end
%Reformat if the column/row ratio is too large: we want a roughly
%square design 
while p(2)/p(1)>2.5
    N=n+1;
    [p,n]=numSubplots(N); % Recursive!
end

end

% End.
