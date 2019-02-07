function info = isquare(A)

%@info:
%! @deftypefn {Function File} {@var{info} =} isquare (@var{A})
%! @anchor{isquare}
%! @sp 1
%! Tests if @var{A} is square matrix.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item A
%! Matlab's array.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item info
%! Integer scalar equal to 1 if @var{A} is a square matrix, 0 otherwise.
%! @end table
%! @end deftypefn
%@eod:

% Copyright (C) 2013-2014 Dynare Team
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

info = 0;
if isequal(ndims(A),2) && isequal(size(A,1),size(A,2))
    info = 1;
end