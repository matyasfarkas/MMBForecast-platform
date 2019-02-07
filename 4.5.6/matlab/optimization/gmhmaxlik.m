function [PostMode, HessianMatrix, Scale, ModeValue] = gmhmaxlik(fun, xinit, Hinit, iscale, bounds, priorstd, gmhmaxlikOptions, OptimizationOptions, varargin)

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

% Set default options

if ~isempty(Hinit)
    gmhmaxlikOptions.varinit = 'previous';
else
    gmhmaxlikOptions.varinit = 'prior';
end

if ~isempty(OptimizationOptions)
    DynareOptionslist = read_key_value_string(OptimizationOptions);
    for i=1:rows(DynareOptionslist)
        switch DynareOptionslist{i,1}
          case 'NumberOfMh'
            gmhmaxlikOptions.iterations = DynareOptionslist{i,2};
          case 'ncov-mh'
            gmhmaxlikOptions.number = DynareOptionslist{i,2};
          case 'nscale-mh'
            gmhmaxlikOptions.nscale = DynareOptionslist{i,2};
          case 'nclimb-mh'
            gmhmaxlikOptions.nclimb = DynareOptionslist{i,2};
          case 'InitialCovarianceMatrix'
            switch DynareOptionslist{i,2}
              case 'previous'
                if isempty(Hinit)
                    error('gmhmaxlik: No previous estimate of the Hessian matrix available! You cannot use the InitialCovarianceMatrix option!')
                else
                    gmhmaxlikOptions.varinit = 'previous';
                end
              case {'prior', 'identity'}
                gmhmaxlikOptions.varinit = DynareOptionslist{i,2};
              otherwise
                error('gmhmaxlik: Unknown value for option ''InitialCovarianceMatrix''!')
            end
          case 'AcceptanceRateTarget'
            gmhmaxlikOptions.target = DynareOptionslist{i,2};
            if gmhmaxlikOptions.target>1 || gmhmaxlikOptions.target<eps
                error('gmhmaxlik: The value of option AcceptanceRateTarget should be a double between 0 and 1!')
            end
          otherwise
            warning(['gmhmaxlik: Unknown option (' DynareOptionslist{i,1}  ')!'])
        end
    end
end

% Evaluate the objective function.
OldModeValue = feval(fun,xinit,varargin{:});

if ~exist('MeanPar','var')
    MeanPar = xinit;
end

switch gmhmaxlikOptions.varinit
  case 'previous'
    CovJump = inv(Hinit);
  case 'prior'
    % The covariance matrix is initialized with the prior
    % covariance (a diagonal matrix) %%Except for infinite variances ;-)
    stdev = priorstd;
    indx = find(isinf(stdev));
    stdev(indx) = ones(length(indx),1)*sqrt(10);
    vars = stdev.^2;
    CovJump = diag(vars);
  case 'identity'
    vars = ones(length(priorstd),1)*0.1;
    CovJump = diag(vars);
  otherwise
    error('gmhmaxlik: This is a bug! Please contact the developers.')
end

OldPostVariance = CovJump;
OldPostMean = xinit;
OldPostMode = xinit;
Scale = iscale;

for i=1:gmhmaxlikOptions.iterations
    if i<gmhmaxlikOptions.iterations
        flag = '';
    else
        flag = 'LastCall';
    end
    [PostMode, PostVariance, Scale, PostMean] = gmhmaxlik_core(fun, OldPostMode, bounds, gmhmaxlikOptions, Scale, flag, MeanPar, OldPostVariance, varargin{:});
    ModeValue = feval(fun, PostMode, varargin{:});
    dVariance = max(max(abs(PostVariance-OldPostVariance)));
    dMean = max(abs(PostMean-OldPostMean));
    skipline()
    printline(58,'=')
    disp(['   Change in the posterior covariance matrix = ' num2str(dVariance) '.'])
    disp(['   Change in the posterior mean = ' num2str(dMean) '.'])
    disp(['   Mode improvement = ' num2str(abs(OldModeValue-ModeValue))])
    disp(['   New value of jscale = ' num2str(Scale)])
    printline(58,'=')
    OldModeValue = ModeValue;
    OldPostMean = PostMean;
    OldPostVariance = PostVariance;
end

HessianMatrix = inv(PostVariance);

skipline()
disp(['Optimal value of the scale parameter = ' num2str(Scale)])
skipline()