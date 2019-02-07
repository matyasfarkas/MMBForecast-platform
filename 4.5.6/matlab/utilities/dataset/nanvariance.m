function variances = nanvariance(data)
% Compute the standard deviation for each observed variable (possibly with missing observations).

%@info:
%! @deftypefn {Function File} {@var{variances} =} nanvariance(@var{data})
%! @anchor{nanvariance}
%! This function computes the variances of the observed variables (possibly with missing observations).
%!
%! @strong{Inputs}
%! @table @var
%! @item datas
%! A T*N array of real numbers.
%! @end table
%!
%! @strong{Outputs}
%! @table @var
%! @item variances
%! A N*1 vector of real numbers
%! @end table
%!
%! @strong{This function is called by:}
%! @ref{descriptive_statistics}.
%!
%! @strong{This function calls:}
%! @ref{ndim}, @ref{demean}, @ref{nandemean}.
%!
%! @strong{Remark 1.} On exit, a new field is appended to the structure: @code{dataset_.descriptive.stdv} is a
%! @tex{n\times 1} vector (where @tex{n} is the number of observed variables as defined by @code{dataset_.info.nvobs}).
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

if isanynan(data)
    variances = transpose(nanmean(bsxfun(@power,nandemean(data),2)));
else
    variances = transpose(mean(bsxfun(@power,demean(data),2)));
end