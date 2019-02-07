function autocov = nanautocovariance(data,order)
% Computes the (multivariate) auto-covariance function of the sample (possibly with missing observations).

%@info:
%! @deftypefn {Function File} {@var{dataset_} =} compute_corr(@var{dataset_},@var{nlag})
%! @anchor{compute_acov}
%! This function computes the (multivariate) auto-covariance function of the sample (possibly with missing observations).
%!
%! @strong{Inputs}
%! @table @var
%! @item data
%! T*N array of real numbers.
%! @item order
%! Integer scalar. The maximum number of lags of the autocovariance function.
%! @end table
%!
%! @strong{Outputs}
%! @table @var
%! @item autocov
%! A N*N*order array of real numbers.
%! @end table
%!
%! @strong{This function is called by:}
%! @ref{descriptive_statistics}.
%!
%! @strong{This function calls:}
%! @ref{ndim}, @ref{nancovariance}, @ref{demean}, @ref{nandemean}.
%!
%! @strong{Remark 1.} On exit, a new field is appended to the structure: @code{dataset_.descriptive.acov} is a
%! @tex{n\times n\times p} array (where @tex{n} is the number of observed variables as defined by @code{dataset_.info.nvobs},
%! and @tex{n} is the maximum number of lags given by the second input @code{nlag}).
%!
%! @strong{Remark 2.} If @code{dataset_.descriptive.cova} does not exist, the covariance matrix is computed prior to the
%! computation of the auto-covariance function.
%!
%! @end deftypefn
%@eod:

% Copyright (C) 2011-2014 Dynare Team
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

n = size(data,2);
missing = isanynan(data);

autocov = zeros(n, n, order);

for lag=1:order
    if missing
        data = nandemean(data);
    else
        data = demean(data);
    end
    for i=1:n
        for j=1:n
            if missing
                autocov(i,j,lag) = nanmean(data((lag+1):end,i).*data(1:end-lag,j));
            else
                autocov(i,j,lag) = mean(data((lag+1):end,i).*data(1:end-lag,j));
            end
        end
    end
end