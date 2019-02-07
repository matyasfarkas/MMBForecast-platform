function c = nandemean(x)
% Removes the mean of each column of a matrix with some NaNs.

%@info:
%! @deftypefn {Function File} {@var{c} =} nandemean (@var{x})
%! @anchor{nandemean}
%! @sp 1
%! This function removes the mean of each column of a matrix with some NaNs.
%! @sp 2
%! @strong{Inputs}
%! @table @var
%! @item x
%! Matlab matrix (T-by-N).
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @table @var
%! @item c
%! Matlab matrix (T-by-N). The demeaned x matrix.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 1
%! @ref{compute_cova}, @ref{compute_acov}, @ref{compute_std}
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
    c = x-nanmean(x);
  case 2
    c = bsxfun(@minus,x,nanmean(x));
  otherwise
    error('descriptive_statistics::nandemean:: This function is not implemented for arrays with dimension greater than two!')
end