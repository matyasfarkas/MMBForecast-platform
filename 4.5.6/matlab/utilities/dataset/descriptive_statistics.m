function dataset_ = descriptive_statistics(dataset_,statistic,varagin)
% Computes various descriptive statistics for the sample and stores them in the structure dataset_.

%@info:
%! @deftypefn {Function File} {@var{dataset_} =} descriptive_statistics(@var{dataset_},@var{statistic})
%! @deftypefn {Function File} {@var{dataset_} =} descriptive_statistics(@var{dataset_},@var{statistic},nlags)
%! @anchor{compute_corr}
%! This function computes various descriptive statistics on the sample (possibly with missing observations).
%!
%! @strong{Inputs}
%! @table @var
%! @item dataset_
%! Dynare structure describing the dataset, built by @ref{initialize_dataset}
%! @item statistic
%! String. The name of the statistic to be computed. Admissible values are:
%!   @table @var
%!   @item 'stdv'
%!   Computes the standard deviation of each observed variable.
%!   @item 'cova'
%!   Computes the covariance matrix of the sample.
%!   @item 'corr'
%!   Computes the correlation matrix of the sample.
%!   @item 'acov'
%!   Computes the (multivariate) auto-covariance function of the sample. In this case a third argument (@code{nlags}) defining the
%!   maximum number of lags is mandatory.
%!   @end table
%! @item nlags
%! Integer scalar. The maximum number of lags when computing the autocovariance function.
%! @end table
%!
%! @strong{Outputs}
%! @table @var
%! @item dataset_
%! Dynare structure describing the dataset, built by @ref{initialize_dataset}
%! @end table
%!
%! @strong{This function is called by:}
%! none.
%!
%! @strong{This function calls:}
%! @ref{compute_stdv}, @ref{compute_cova}, @ref{compute_corr}, @ref{compute_acov}.
%!
%! @strong{Remark 1.} On exit, a new field containing the computed statistics is appended to the structure.
%!
%! @end deftypefn
%@eod:

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

% Original author: stephane DOT adjemian AT univ DASH lemans DOT fr


if strcmpi(statistic,'stdv')
    dataset_ = compute_std(dataset_)
end

if strcmpi(statistic,'cova')
    dataset_ = compute_cova(dataset_)
end

if strcmpi(statistic,'corr')
    dataset_ = compute_cova(dataset_)
end

if strcmpi(statistic,'acov')
    if nargin==2
        nlag = 10;
    else
        nlag = varargin{1};
    end
    dataset_ = compute_acov(dataset_,nlag);
end