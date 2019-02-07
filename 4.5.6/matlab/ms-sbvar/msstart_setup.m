%function []= msstart_setup(options_)

% Copyright (C) 2011-2017 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

% ** ONLY UNDER UNIX SYSTEM
%path(path,'/usr2/f1taz14/mymatlab')



%===========================================
% Exordium I
%===========================================
format short g     % format
                   %
                   %options_.ms.freq = 4;   % quarters or months
                   %options_.ms.initial_year=1959;   % beginning of the year
                   %options_.ms.initial_subperiod=1;    % begining of the quarter or month
                   %options_.ms.final_year=2005;   % final year
                   %options_.ms.final_subperiod=4;    % final month or quarter
nData=(options_.ms.final_year-options_.ms.initial_year)*options_.ms.freq + (options_.ms.final_subperiod-options_.ms.initial_subperiod+1);
% total number of the available data -- this is all you have

%*** Load data and series
%load datainf_argen.prn      % the default name for the variable is "options_.ms.data".
%load datacbogdpffr.prn
%options_.ms.data = datacbogdpffr;
%clear datacbogdpffr;
[nt,ndv]=size(options_.data);
if nt < nData
    error('The declared sample is longer than the available data')
end
%--------
%1    CBO output gap --  log(x_t)-log(x_t potential)
%2    GDP deflator -- (P_t/P_{t-1})^4-1.0
%2    FFR/100.
options_.ms.vlist = [1:length(options_.varobs)];    % 1: U; 4: PCE inflation.
options_.ms.varlist=cellstr(options_.varobs');
options_.ms.log_var = sort(varlist_indices(options_.ms.vlistlog,char(options_.varobs)));   % subset of "options_.ms.vlist.  Variables in log level so that differences are in **monthly** growth, unlike R and U which are in annual percent (divided by 100 already).
options_.ms.percent_var =setdiff(options_.ms.vlist,options_.ms.log_var);
%options_.ms.restriction_fname='ftd_upperchol3v';   %Only used by msstart2.m.
ylab = options_.ms.varlist;
xlab = options_.ms.varlist;

%----------------
nvar = length(options_.varobs);   % number of endogenous variables
nlogeno = length(options_.ms.log_var);  % number of endogenous variables in options_.ms.log_var
npereno = length(options_.ms.percent_var);  % number of endogenous variables in options_.ms.percent_var
if (nvar~=(nlogeno+npereno))
    skipline()
    warning('Check xlab, nlogeno or npereno to make sure of endogenous variables in options_.ms.vlist')
    disp('Press ctrl-c to abort')
    return
elseif (nvar==length(options_.ms.vlist))
    nexo=1;    % only constants as an exogenous variable.  The default setting.
elseif (nvar<length(options_.ms.vlist))
    nexo=length(options_.ms.vlist)-nvar+1;
else
    skipline()
    warning('Make sure there are only nvar endogenous variables in options_.ms.vlist')
    disp('Press ctrl-c to abort')
    return
end


%------- A specific sample is considered for estimation -------
yrStart=options_.ms.initial_year;
qmStart=options_.ms.initial_subperiod;
yrEnd=options_.ms.final_year;
qmEnd=options_.ms.final_subperiod;
%options_.forecast = 4;   % number of years for forecasting
if options_.forecast<1
    error('To be safe, the number of forecast years should be at least 1')
end
forelabel = [num2str(yrEnd) ':' num2str(qmEnd) ' Forecast'];

nSample=(yrEnd-yrStart)*options_.ms.freq + (qmEnd-qmStart+1);
if qmEnd==options_.ms.freq
    E1yrqm = [yrEnd+1 1];  % first year and quarter (month) after the sample
else
    E1yrqm = [yrEnd qmEnd+1];  % first year and quarter (month) after the sample
end
E2yrqm = [yrEnd+options_.forecast qmEnd];   % end at the last month (quarter) of a calendar year after the sample
[fdates,nfqm]=fn_calyrqm(options_.ms.freq,E1yrqm,E2yrqm);   % forecast dates and number of forecast dates
[sdates,nsqm] = fn_calyrqm(options_.ms.freq,[yrStart qmStart],[yrEnd qmEnd]);
% sdates: dates for the whole sample (including options_.ms.nlags)
if nSample~=nsqm
    warning('Make sure that nSample is consistent with the size of sdates')
    disp('Hit any key to continue, or ctrl-c to abort')
    pause
end
imstp = 4*options_.ms.freq;    % <<>>  impulse responses (4 years)
nayr = 4; %options_.forecast;  % number of years before forecasting for plotting.


%------- Prior, etc. -------
%options_.ms.nlags = 4;        % number of options_.ms.nlags
%options_.ms.cross_restrictions = 0;   % 1: cross-A0-and-A+ restrictions; 0: options_.ms.restriction_fname is all we have
% Example for indxOres==1: restrictions of the form P(t) = P(t-1).
%options_.ms.contemp_reduced_form = 0;  % 1: contemporaneous recursive reduced form; 0: restricted (non-recursive) form
%options_.ms.real_pseudo_forecast = 0;  % 1: options_.ms.real_pseudo_forecast forecasts; 0: real time forecasts
%options_.ms.bayesian_prior = 1;  % 1: Bayesian prior; 0: no prior
indxDummy = options_.ms.bayesian_prior;  % 1: add dummy observations to the data; 0: no dummy added.
                                         %options_.ms.dummy_obs = 0;  % No dummy observations for xtx, phi, fss, xdatae, etc.  Dummy observations are used as an explicit prior in fn_rnrprior_covres_dobs.m.
                                         %if indxDummy
                                         %   options_.ms.dummy_obs=nvar+1;         % number of dummy observations
                                         %else
                                         %   options_.ms.dummy_obs=0;    % no dummy observations
                                         %end
                                         %=== The following mu is effective only if options_.ms.bayesian_prior==1.

mu = options_.ms.coefficients_prior_hyperparameters;

%   mu(1): overall tightness and also for A0;
%   mu(2): relative tightness for A+;
%   mu(3): relative tightness for the constant term;
%   mu(4): tightness on lag decay;  (1)
%   mu(5): weight on nvar sums of coeffs dummy observations (unit roots);
%   mu(6): weight on single dummy initial observation including constant
%           (cointegration, unit roots, and stationarity);
%
%
hpmsmd = [0.0; 0.0];
indxmsmdeqn = [0; 0; 0; 0];  %This option disenable using this in fn_rnrprior_covres_dobs.m


tdf = 3;          % degrees of freedom for t-dist for initial draw of the MC loop
nbuffer = 1000;        % a block or buffer of draws (buffer) that is saved to the disk (not memory)
ndraws1=1*nbuffer;         % 1st part of Monte Carlo draws
ndraws2=10*ndraws1;        % 2nd part of Monte Carlo draws
                           % seednumber = options_.DynareRandomStreams.seed; %7910;    %472534;   % if 0, random state at each clock time
                           %            % good one 420 for [29 45], [29 54]
                           % if seednumber
                           %    randn('state',seednumber);
                           %    rand('state',seednumber);
                           % else
                           %    randn('state',fix(100*sum(clock)));
                           %    rand('state',fix(100*sum(clock)));
                           % end
                           %  nstarts=1         % number of starting points
                           %  imndraws = nstarts*ndraws2;   % total draws for impulse responses or forecasts
                           %<<<<<<<<<<<<<<<<<<<
