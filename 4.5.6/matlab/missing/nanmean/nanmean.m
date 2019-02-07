function y = nanmean(x)
% Computes the mean of each column of a matrix with some NaNs.

%@info:
%! @deftypefn {Function File} {@var{y} =} nanmean (@var{x})
%! @anchor{nanmean}
%! @sp 1
%! Computes the mean of each column of a matrix with some NaNs.
%! @sp 2
%! @strong{Inputs}
%! @table @ @var
%! @item x
%! Matlab matrix (T-by-N).
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @table @ @var
%! @item y
%! Matlab vector (1-by-N), the mean.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 1
%! @ref{compute_cova}, @ref{compute_acov}, @ref{compute_std}, @ref{nandemean}
%! @sp 2
%! @strong{This function calls:}
%! @sp 1
%! @ref{ndim}
%!
%! @end deftypefn
%@eod:

% Copyright (C) 2011 Dynare Team
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

% AUTHOR(S) stephane DOT adjemian AT univ DASH lemans DOT fr

switch ndim(x)
  case 1
    y = mean(x(find(~isnan(x))));
  case 2
    y = NaN(1,size(x,2));
    for i = 1:size(x,2)
        y(i) = mean(x(find(~isnan(x(:,i))),i));
    end
  otherwise
    error('descriptive_statistics::nanmean:: This function is not implemented for arrays with dimension greater than two!')
end