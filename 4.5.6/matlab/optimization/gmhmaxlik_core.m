function [PostMod,PostVar,Scale,PostMean] = gmhmaxlik_core(ObjFun,xparam1,mh_bounds,options,iScale,info,MeanPar,VarCov,varargin)

% (Dirty) Global minimization routine of (minus) a likelihood (or posterior density) function.
%
% INPUTS
%   o ObjFun     [char]     string specifying the name of the objective function.
%   o xparam1    [double]   (p*1) vector of parameters to be estimated.
%   o mh_bounds  [double]   (p*2) matrix defining lower and upper bounds for the parameters.
%   o options    [structure] options for the optimization algorithm (options_.gmhmaxlik).
%   o iScale     [double]   scalar specifying the initial of the jumping distribution's scale parameter.
%   o info       [char]     string, empty or equal to 'LastCall'.
%   o MeanPar    [double]   (p*1) vector specifying the initial posterior mean.
%   o VarCov     [double]   (p*p) matrix specifying the initial posterior covariance matrix.
%   o gend       [integer]  scalar specifying the number of observations ==> varargin{1}.
%   o data       [double]   (T*n) matrix of data ==> varargin{2}.
%
% OUTPUTS
%   o PostMod    [double]   (p*1) vector, evaluation of the posterior mode.
%   o PostVar    [double]   (p*p) matrix, evaluation of the posterior covariance matrix.
%   o Scale      [double]   scalar specifying the scale parameter that should be used in
%                           an eventual metropolis-hastings algorithm.
%   o PostMean   [double]   (p*1) vector, evaluation of the posterior mean.
%
% ALGORITHM
%   Metropolis-Hastings with an constantly updated covariance matrix for
%   the jump distribution. The posterior mean, variance and mode are
%   updated (in step 2) with the following rules:
%
%   \[
%       \mu_t = \mu_{t-1} + \frac{1}{t}\left(\theta_t-\mu_{t-1}\right)
%   \]
%
%   \[
%       \Sigma_t = \Sigma_{t-1} + \mu_{t-1}\mu_{t-1}'-\mu_{t}\mu_{t}' +
%                  \frac{1}{t}\left(\theta_t\theta_t'-\Sigma_{t-1}-\mu_{t-1}\mu_{t-1}'\right)
%   \]
%
%   and
%
%   \[
%       \mathrm{mode}_t = \left\{
%                       \begin{array}{ll}
%                         \theta_t, & \hbox{if } p(\theta_t|\mathcal Y) > p(\mathrm{mode}_{t-1}|\mathcal Y) \\
%                         \mathrm{mode}_{t-1}, & \hbox{otherwise.}
%                       \end{array}
%                     \right.
%   \]
%
%   where $t$ is the iteration, $\mu_t$ the estimate of the posterior mean
%   after $t$ iterations, $\Sigma_t$ the estimate of the posterior
%   covariance matrix after $t$ iterations, $\mathrm{mode}_t$ is the
%   evaluation of the posterior mode after $t$ iterations and
%   $p(\theta_t|\mathcal Y)$ is the posterior density of parameters
%   (specified by the user supplied function "fun").
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright (C) 2006-2017 Dynare Team
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

global bayestopt_ estim_params_ options_

options_.lik_algo = 1;
npar = length(xparam1);

NumberOfIterations = options.number;
MaxNumberOfTuningSimulations   = options.nscale;
MaxNumberOfClimbingSimulations = options.nclimb;
AcceptanceTarget               = options.target;

CovJump = VarCov;
ModePar = xparam1;

%% [1] I tune the scale parameter.
hh = dyn_waitbar(0,'Tuning of the scale parameter...');
set(hh,'Name','Tuning of the scale parameter.');
j = 1; jj  = 1;
isux = 0; jsux = 0; test = 0;
ix2 = ModePar;% initial condition!
ilogpo2 = - feval(ObjFun,ix2,varargin{:});% initial posterior density
mlogpo2 = ilogpo2;
try
    dd = transpose(chol(CovJump));
catch
    dd = eye(length(CovJump));
end
while j<=MaxNumberOfTuningSimulations
    proposal = iScale*dd*randn(npar,1) + ix2;
    if all(proposal > mh_bounds(:,1)) && all(proposal < mh_bounds(:,2))
        logpo2 = - feval(ObjFun,proposal,varargin{:});
    else
        logpo2 = -inf;
    end
    % I move if the proposal is enough likely...
    if logpo2 > -inf && log(rand) < logpo2 - ilogpo2
        ix2 = proposal;
        if logpo2 > mlogpo2
            ModePar = proposal;
            mlogpo2 = logpo2;
        end
        ilogpo2 = logpo2;
        isux = isux + 1;
        jsux = jsux + 1;
    end% ... otherwise I don't move.
    prtfrc = j/MaxNumberOfTuningSimulations;
    if mod(j, 10)==0
        dyn_waitbar(prtfrc,hh,sprintf('Acceptance ratio [during last 500]: %f [%f]',isux/j,jsux/jj));
    end
    if  j/500 == round(j/500)
        test1 = jsux/jj;
        cfactor = test1/AcceptanceTarget;
        if cfactor>0
            iScale = iScale*cfactor;
        else
            iScale = iScale/10;
        end
        jsux = 0; jj = 0;
        if cfactor>0.90 && cfactor<1.10
            test = test+1;
        end
        if test>4
            break
        end
    end
    j = j+1;
    jj = jj + 1;
end

dyn_waitbar_close(hh);
%% [2] One block metropolis, I update the covariance matrix of the jumping distribution
hh = dyn_waitbar(0,'Metropolis-Hastings...');
set(hh,'Name','Estimation of the posterior covariance...'),
j = 1;
isux = 0;
ilogpo2 = - feval(ObjFun,ix2,varargin{:});
while j<= NumberOfIterations
    proposal = iScale*dd*randn(npar,1) + ix2;
    if all(proposal > mh_bounds(:,1)) && all(proposal < mh_bounds(:,2))
        logpo2 = - feval(ObjFun,proposal,varargin{:});
    else
        logpo2 = -inf;
    end
    % I move if the proposal is enough likely...
    if logpo2 > -inf && log(rand) < logpo2 - ilogpo2
        ix2 = proposal;
        if logpo2 > mlogpo2
            ModePar = proposal;
            mlogpo2 = logpo2;
        end
        ilogpo2 = logpo2;
        isux = isux + 1;
        jsux = jsux + 1;
    end% ... otherwise I don't move.
    prtfrc = j/NumberOfIterations;
    if mod(j, 10)==0
        dyn_waitbar(prtfrc,hh,sprintf('Acceptance ratio: %f',isux/j));
    end
    % I update the covariance matrix and the mean:
    oldMeanPar = MeanPar;
    MeanPar = oldMeanPar + (1/j)*(ix2-oldMeanPar);
    CovJump = CovJump + oldMeanPar*oldMeanPar' - MeanPar*MeanPar' + ...
              (1/j)*(ix2*ix2' - CovJump - oldMeanPar*oldMeanPar');
    j = j+1;
end
dyn_waitbar_close(hh);
PostVar = CovJump;
PostMean = MeanPar;
%% [3 & 4] I tune the scale parameter (with the new covariance matrix) if
%% this is the last call to the routine, and I climb the hill (without
%% updating the covariance matrix)...
if strcmpi(info,'LastCall')

    hh = dyn_waitbar(0,'Tuning of the scale parameter...');
    set(hh,'Name','Tuning of the scale parameter.'),
    j = 1; jj  = 1;
    isux = 0; jsux = 0;
    test = 0;
    ilogpo2 = - feval(ObjFun,ix2,varargin{:});% initial posterior density
    dd = transpose(chol(CovJump));
    while j<=MaxNumberOfTuningSimulations
        proposal = iScale*dd*randn(npar,1) + ix2;
        if all(proposal > mh_bounds(:,1)) && all(proposal < mh_bounds(:,2))
            logpo2 = - feval(ObjFun,proposal,varargin{:});
        else
            logpo2 = -inf;
        end
        % I move if the proposal is enough likely...
        if logpo2 > -inf && log(rand) < logpo2 - ilogpo2
            ix2 = proposal;
            if logpo2 > mlogpo2
                ModePar = proposal;
                mlogpo2 = logpo2;
            end
            ilogpo2 = logpo2;
            isux = isux + 1;
            jsux = jsux + 1;
        end% ... otherwise I don't move.
        prtfrc = j/MaxNumberOfTuningSimulations;
        if mod(j, 10)==0
            dyn_waitbar(prtfrc,hh,sprintf('Acceptance ratio [during last 1000]: %f [%f]',isux/j,jsux/jj));
        end
        if j/1000 == round(j/1000)
            test1 = jsux/jj;
            cfactor = test1/AcceptanceTarget;
            iScale = iScale*cfactor;
            jsux = 0; jj = 0;
            if cfactor>0.90 && cfactor<1.10
                test = test+1;
            end
            if test>4
                break
            end
        end
        j = j+1;
        jj = jj + 1;
    end
    dyn_waitbar_close(hh);
    Scale = iScale;
    %%
    %% Now I climb the hill
    %%
    if options.nclimb
        hh = dyn_waitbar(0,' ');
        set(hh,'Name','Now I am climbing the hill...'),
        j = 1; jj  = 1;
        jsux = 0;
        test = 0;
        while j<=MaxNumberOfClimbingSimulations
            proposal = iScale*dd*randn(npar,1) + ModePar;
            if all(proposal > mh_bounds(:,1)) && all(proposal < mh_bounds(:,2))
                logpo2 = - feval(ObjFun,proposal,varargin{:});
            else
                logpo2 = -inf;
            end
            if logpo2 > mlogpo2% I move if the proposal is higher...
                ModePar = proposal;
                mlogpo2 = logpo2;
                jsux = jsux + 1;
            end% otherwise I don't move...
            prtfrc = j/MaxNumberOfClimbingSimulations;
            if mod(j, 10)==0
                dyn_waitbar(prtfrc,hh,sprintf('%f Jumps / MaxStepSize %f',jsux,sqrt(max(diag(iScale*CovJump)))));
            end
            if  j/200 == round(j/200)
                if jsux<=1
                    test = test+1;
                else
                    test = 0;
                end
                jsux = 0;
                jj = 0;
                if test>4% If I do not progress enough I reduce the scale parameter
                         % of the jumping distribution (cooling down the system).
                    iScale = iScale/1.10;
                end
                if sqrt(max(diag(iScale*CovJump)))<10^(-4)
                    break% Steps are too small!
                end
            end
            j = j+1;
            jj = jj + 1;
        end
        dyn_waitbar_close(hh);
    end%climb
else
    Scale = iScale;
end
PostMod = ModePar;