function ms_sbvar_setup(options_)
% function ms_sbvar_setup(options_)
% does the general file initialization for ms sbvar
%
% INPUTS
%    options_:    (struct)    options
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2003-2017 Dynare Team
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

options_.data = read_variables(options_.datafile, ...
                               options_.varobs, [], options_.xls_sheet, options_.xls_range);
[options_.ms.final_year,options_.ms.final_subperiod] = check_datafile_years_assigned(options_);

if options_.ms.upper_cholesky
    if options_.ms.lower_cholesky
        error(['Upper Cholesky and lower Cholesky decomposition can''t be ' ...
               'requested at the same time!'])
    else
        options_.ms.restriction_fname = 'upper_cholesky';
    end
elseif options_.ms.lower_cholesky
    options_.ms.restriction_fname = 'lower_cholesky';
elseif ~isempty(options_.ms.Qi) && ~isempty(options_.ms.Ri)
    options_.ms.restriction_fname = 'exclusions';
else
    options_.ms.restriction_fname = 0;
end

%==========================================================================
%== Markov Process Specification File
%==========================================================================
markov_file = [options_.ms.output_file_tag '_markov_file.dat'];

%==========================================================================
%== BVAR prior
%==========================================================================

%=== The following mu is effective only if indxPrior==1.
%mu = zeros(6,1);   % hyperparameters
if length(options_.ms.coefficients_prior_hyperparameters) ~= 6
    error('When specifying the coefficients_prior_hyperparameters, you must pass a vector of 6 numbers')
end
mu = options_.ms.coefficients_prior_hyperparameters;
mu = reshape(mu,1,6);

% mu(1): overall tightness for A0 and Aplus
% mu(2): relative tightness for Aplus
% mu(3): relative tightness for the constant term
% mu(4): tightness on lag decay.  (1.2 - 1.5 faster decay produces better
% inflation forecasts
% mu(5): weight on nvar sums of coeffs dummy observations (unit roots).
% mu(6): weight on single dummy initial observation including constant
% (cointegration, unit roots, and stationarity).

% Alpha on p. 66 for squared time-varying structural shock lambda.
galp = options_.ms.alpha;

% Beta on p. 66 for squared time-varying structural shock lambda.
gbeta = options_.ms.beta;

% Case 3 (no state change across options_.ms.nlags (l) but allows all variables for a
% given lag to switch states). Normal prior variance for glamda
% (nvar-by-nvar for each state) for different variables in lagged D+.  See
% p.71v.
gsig2_lmdm = options_.ms.gsig2_lmdm;


%==========================================================================
%== Data
%==========================================================================
% Read in data to produce rectangular array named xdd.  Each column is one
% data series.
xdd=options_.data;

% Information about timing of the data for consistancy checks
% quarters (4) or months (12)
q_m = options_.ms.freq;
% beginning year in data set
yrBin=options_.ms.initial_year;
% beginning quarter or month in data set
%options_.ms.initial_subperiod = 1;
qmBin=options_.ms.initial_subperiod;
% final year in data set
yrFin=options_.ms.final_year;
% final month or quarter in data set
qmFin=options_.ms.final_subperiod;
% first year to use in estimation
yrStart=options_.ms.initial_year;
% first quarter or month to use in estimation
qmStart=options_.ms.initial_subperiod;
% last year to use in estimation
yrEnd=options_.ms.final_year;
% last quater or month to use in estimation
qmEnd=options_.ms.final_subperiod;
% Log variables in xdd
logindx = [];

% Convert percent to decimal in xdd
pctindx = [];

% Select the variable to use and rearrange columns if desired
%vlist = [3 1 2];
%options_.ms.vlist = [1 2 3];
options_.ms.vlist = 1:length(options_.varobs);
vlist1=options_.ms.vlist;

%==========================================================================
%==========================================================================
%==========================================================================
%== Beginning of code.  Modify below at own risk.
%==========================================================================

% options that may at some point become user specified
%indxC0Pres = 0;    % 1: cross-A0-and-A+ restrictions; 0: idfile_const is all we have
indxC0Pres =options_.ms.cross_restrictions;
% Example for indxOres==1: restrictions of the form P(t) = P(t-1).
%Rform = 0;         % 1: contemporaneous recursive reduced form; 0: restricted (non-recursive) form
Rform =options_.ms.contemp_reduced_form;
% % % Pseudo = 0;  % 1: Pseudo forecasts; 0: real time forecasts
%indxPrior = 1;     % 1: Bayesian prior; 0: no prior
indxPrior =options_.ms.bayesian_prior;
%indxDummy = indxPrior;  % 1: add dummy observations to the data; 0: no dummy added.
indxDummy = options_.ms.bayesian_prior;
%ndobs = 0;         % No dummy observations for xtx, phi, fss, xdatae, etc.  Dummy observations are used as an explicit prior in fn_rnrprior_covres_dobs.m.
%ndobs =options_.ms.dummy_obs;
%if indxDummy
%   ndobs=nvar+1;         % number of dummy observations
%else
%   ndobs=0;    % no dummy observations
%end

%
hpmsmd = [0.0; 0.0];
indxmsmdeqn = [0; 0; 0; 0];  %This option disenable using this in fn_rnrprior_covres_dobs.m

nStates = -1;

%==========================================================================
%== Create initialization file
%==========================================================================

%======================================================================
%== Check and setup data
%======================================================================
% log data
xdd(:,logindx) = log(xdd(:,logindx));

% convert percentage to decimal
xdd(:,pctindx)=.01*xdd(:,pctindx);

if (q_m ~= 12) && (q_m ~= 4)
    disp('Warning: data must be monthly or quarterly!')
    return
end

% number of data points
nData=(yrFin-yrBin)*q_m + (qmFin-qmBin+1);
% number of data points in estimation sample
nSample=(yrEnd-yrStart)*q_m + (qmEnd-qmEnd+1);
% number of periods not used at beginning of data (non-negative number)
nStart=(yrStart-yrBin)*q_m + (qmStart-qmBin);
% number of periods not used at end of data (non-positive number)
nEnd=(yrEnd-yrFin)*q_m + (qmEnd-qmFin);

if (nEnd > 0) || (nStart < 0)
    disp('Warning: desired estimation period not in data set!')
    return
end
if (nSample <= 0)
    disp('Warning: no data points in estimation period!')
    return
end

% reorder variables and create estimation data set
xdgel=xdd(nStart+1:nData+nEnd,vlist1);

% bad data points
baddata = find(isnan(xdgel));
if ~isempty(baddata)
    disp('Warning: some data for estimation period are unavailable.')
    return
end

% set nvar and nexo
nvar=size(xdgel,2);
nexo=1;

% Arranged data information, WITHOUT dummy obs when 0 after mu is used.
% See fn_rnrprior_covres_dobs.m for using the dummy observations as part of
% an explicit prior.
[xtx,xty,yty,fss,phi,y,ncoef,xr,Bh] = fn_dataxy(nvar,options_.ms.nlags,xdgel,mu,0,nexo);


%======================================================================
%== Linear Restrictions
%======================================================================
if Rform
    Ui=cell(nvar,1); Vi=cell(ncoef,1);
    for kj=1:nvar
        Ui{kj} = eye(nvar);  Vi{kj} = eye(ncoef);
    end
else
    [Ui,Vi,n0,np,ixmC0Pres] = feval(options_.ms.restriction_fname,nvar,nexo,options_.ms);
    if min(n0)==0
        disp('A0: restrictions give no free parameters in one of equations')
        return
    elseif min(np)==0
        disp('Aplus: Restrictions in give no free parameters in one of equations')
        return
    end
end


%======================================================================
%== Estimation
%======================================================================
if indxPrior
    %*** Obtains asymmetric prior (with no linear restrictions) with dummy observations as part of an explicit prior (i.e,
    %      reflected in Hpmulti and Hpinvmulti).  See Forecast II, pp.69a-69b for details.
    if 1  % Liquidity effect prior on both MS and MD equations.
        [Pi,H0multi,Hpmulti,H0invmulti,Hpinvmulti] = fn_rnrprior_covres_dobs(nvar,q_m,options_.ms.nlags,xdgel,mu,indxDummy,hpmsmd,indxmsmdeqn);
    else
        [Pi,H0multi,Hpmulti,H0invmulti,Hpinvmulti] = fn_rnrprior(nvar,q_m,options_.ms.nlags,xdgel,mu);
    end

    %*** Combines asymmetric prior with linear restrictions
    [Ptld,H0invtld,Hpinvtld] = fn_rlrprior(Ui,Vi,Pi,H0multi,Hpmulti,nvar);

    %*** Obtains the posterior matrices for estimation and inference
    [Pmat,H0inv,Hpinv] = fn_rlrpostr(xtx,xty,yty,Ptld,H0invtld,Hpinvtld,Ui,Vi);
else
    %*** Obtain the posterior matrices for estimation and inference
    [Pmat,H0inv,Hpinv] = fn_dlrpostr(xtx,xty,yty,Ui,Vi);
end

if Rform
    %*** Obtain the ML estimate
    A0hatinv = chol(H0inv{1}/fss);   % upper triangular but lower triangular choleski
    A0hat=inv(A0hatinv);

    Aphat = Pmat{1}*A0hat;
else
    %*** Obtain the ML estimate
    %   load idenml
    x = 10*rand(sum(n0),1);
    H0 = eye(sum(n0));
    crit = 1.0e-9;
    nit = 10000;
    %
    [fhat,xhat,grad,Hhat,itct,fcount,retcodehat] = csminwel('fn_a0freefun',x,H0,'fn_a0freegrad',crit,nit,Ui,nvar,n0,fss,H0inv);

    A0hat = fn_tran_b2a(xhat,Ui,nvar,n0);

    xhat = fn_tran_a2b(A0hat,Ui,nvar,n0);
    [Aphat,ghat] = fn_gfmean(xhat,Pmat,Vi,nvar,ncoef,n0,np);
    if indxC0Pres
        Fhatur0P = Fhat;  % ur: unrestriced across A0 and A+
        for ki = 1:size(ixmC0Pres,1)   % loop through the number of equations in which
                                       % cross-A0-A+ restrictions occur. See St. Louis Note p.5.
            ixeq = ixmC0Pres{ki}(1,1);   % index for the jth equation in consideration.
            Lit = Vi{ixeq}(ixmC0Pres{ki}(:,2),:);  % transposed restriction matrix Li
                                                   % V_j(i,:) in f_j(i) = V_j(i,:)*g_j
            ci = ixmC0Pres{ki}(:,4) .* A0hat(ixmC0Pres{ki}(:,3),ixeq);
            % s * a_j(h) in the restriction f_j(i) = s * a_j(h).
            LtH = Lit/Hpinv{ixeq};
            HLV = LtH'/(LtH*Lit');
            gihat = Vi{ixeq}'*Fhatur0P(:,ixeq);
            Aphat(:,ixeq) = Vi{ixeq}*(gihat + HLV*(ci-Lit*gihat));
        end
    end
end


%======================================================================
%== Create matlab initialization file
%======================================================================
matlab_filename = ['matlab_',options_.ms.output_file_tag,'.prn'];
fidForC = fopen(matlab_filename,'w');

fprintf(fidForC,'\n%s\n','//== gxia: alpha parameter for gamma prior of xi ==//');
fprintf(fidForC,' %20.15f ', galp);
fprintf(fidForC, '\n\n');

fprintf(fidForC,'\n%s\n','//== gxib: beta parameter for gamma prior of xi ==//');
fprintf(fidForC,' %20.15f ', gbeta);
fprintf(fidForC, '\n\n');

fprintf(fidForC,'\n%s\n','//== glamdasig: sigma parameter for normal prior of lamda ==//');
fprintf(fidForC,' %20.15f ', sqrt(gsig2_lmdm));
fprintf(fidForC, '\n\n');

%=== lags, nvar, nStates, sample size (excluding options_.ms.nlags where, with dummyies, fss=nSample-options_.ms.nlags+ndobs).
fprintf(fidForC,'\n%s\n','//== lags, nvar, nStates, T ==//');
fprintf(fidForC,' %d  %d  %d  %d\n\n\n',options_.ms.nlags, nvar, nStates, fss);

%=== A0hat nvar-by-nvar from the constant VAR.
fprintf(fidForC,'\n%s\n','//== A0hat: nvar-by-nvar ==//');
indxFloat = 1;
xM = A0hat;
nrows = nvar;
ncols = nvar;
fn_fprintmatrix(fidForC, xM, nrows, ncols, indxFloat)

%=== Aphat ncoef-by-nvar from the constant VAR.
%=== Each column of Aphat is in the order of [nvar variables for 1st lag, ..., nvar variables for last lag, constant term].
fprintf(fidForC,'\n%s\n','//== Aphat: ncoef(lags*nvar+1)-by-nvar ==//');
indxFloat = 1;
xM = Aphat;
nrows = ncoef;
ncols = nvar;
fn_fprintmatrix(fidForC, xM, nrows, ncols, indxFloat)

%=== n0const: nvar-by-1, whose ith element represents the number of free A0 parameters in ith equation for the case of constant parameters.
fprintf(fidForC,'\n%s\n','//== n0const: nvar-by-1 ==//');
indxFloat = 0;
xM = n0;
nrows = 1;
ncols = nvar;
fn_fprintmatrix(fidForC, xM', nrows, ncols, indxFloat)

%=== npconst: nvar-by-1, whose ith element represents the number of free A+ parameters in ith equation for the case of constant parameters.
fprintf(fidForC,'\n%s\n','//== npconst: nvar-by-1 ==//');
indxFloat = 0;
xM = np;
nrows = 1;
ncols = nvar;
fn_fprintmatrix(fidForC, xM', nrows, ncols, indxFloat)

%=== Specification
fprintf(fidForC,'\n%s','//== Specification (0=default  1=Sims-Zha  2=Random Walk) ==//');
fprintf(fidForC,'\n%d\n\n',options_.ms.specification);

%=== Uiconst: nvar-by-1 cell.  In each cell, nvar-by-qi orthonormal basis for the null of the ith
%           equation contemporaneous restriction matrix where qi is the number of free parameters.
%           With this transformation, we have ai = Ui*bi or Ui'*ai = bi where ai is a vector
%           of total original parameters and bi is a vector of free parameters. When no
%           restrictions are imposed, we have Ui = I.  There must be at least one free
%           parameter left for the ith equation.
fprintf(fidForC,'\n%s\n','//== Uiconst: cell(nvar,1) and nvar-by-n0const(i) for the ith cell (equation) ==//');
for i_=1:nvar
    fn_fprintmatrix(fidForC, Ui{i_}, nvar, n0(i_), 1);
end

%=== Viconst: nvar-by-1 cell.  In each cell, k-by-ri orthonormal basis for the null of the ith
%           equation lagged restriction matrix where k is a total of exogenous variables and
%           ri is the number of free parameters. With this transformation, we have fi = Vi*gi
%           or Vi'*fi = gi where fi is a vector of total original parameters and gi is a
%           vector of free parameters. There must be at least one free parameter left for
%           the ith equation.
fprintf(fidForC,'\n%s\n','//== Viconst: cell(nvar,1) and ncoef-by-n0const(i) for the ith cell (equation) ==//');
for i_=1:nvar
    fn_fprintmatrix(fidForC, Vi{i_}, ncoef, np(i_), 1);
end

%=== H0barconstcell: cell(nvar,1) (equations) and n-by-n for each cell (equaiton).
%=== H0barconst:  prior covariance matrix for each column of A0 under asymmetric prior (including SZ dummy obs.) with NO linear restrictions imposed yet.
fprintf(fidForC,'\n%s\n','//== H0barconstcell: cell(nvar,1) and n-by-n for the ith cell (equation) ==//');
for i_=1:nvar
    fn_fprintmatrix(fidForC, H0multi(:,:,i_), nvar, nvar, 1);
end

%=== Hpbarconstcell: cell(nvar,1) (equations) and ncoef-by-ncoef for each cell (equaiton).
%=== Hpbarconst:  prior covariance matrix for each column of A+ under asymmetric prior (including SZ dummy obs.) with NO linear restrictions imposed yet.
fprintf(fidForC,'\n%s\n','//== Hpbarconstcell: cell(nvar,1) and ncoef-by-ncoef for the ith cell (equation) ==//');
for i_=1:nvar
    fn_fprintmatrix(fidForC, Hpmulti(:,:,i_), ncoef, ncoef, 1);
end

%=== phi:  X; T-by-k; column: [nvar for 1st lag, ..., nvar for last lag, other exogenous terms, const term]
fprintf(fidForC,'\n%s\n','//== Xright -- X: T-by-ncoef ==//');
xM = phi;
nrows = fss;
ncols = ncoef;
for ki=1:nrows
    for kj=1:ncols
        fprintf(fidForC,' %20.15f ',xM((kj-1)*nrows+ki));
        if (kj==ncols)
            fprintf(fidForC,'\n');
        end
    end
    if (ki==nrows)
        fprintf(fidForC,'\n\n');
    end
end

%=== y:    Y: T-by-nvar where T=fss
fprintf(fidForC,'\n%s\n','//== Yleft -- Y: T-by-nvar ==//');
xM = y;
nrows = fss;
ncols = nvar;
for ki=1:nrows
    for kj=1:ncols
        fprintf(fidForC,' %20.15f ',xM((kj-1)*nrows+ki));
        if (kj==ncols)
            fprintf(fidForC,'\n');
        end
    end
    if (ki==nrows)
        fprintf(fidForC,'\n\n');
    end
end

fclose(fidForC);

%======================================================================
%== Create C initialization filename
%======================================================================
ms_write_markov_file(markov_file,options_)
create_init_file = [matlab_filename,' ',markov_file,' ',options_.ms.file_tag];
[err] = ms_sbvar_create_init_file(create_init_file);
mexErrCheck('ms_sbvar_create_init_file',err);
end
