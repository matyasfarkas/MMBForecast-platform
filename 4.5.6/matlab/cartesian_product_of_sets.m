function cprod = cartesian_product_of_sets(varargin)
% Computes the cartesian product.

%@info:
%! @deftypefn {Function File} {@var{cprod} =} cartesian_product_of_sets (@var{a},@var{b}, ...)
%! @anchor{cartesian_product_of_sets}
%! @sp 1
%! Computes A_1 * A_2 * .... * A_n with a generic set A_i = {e_1,e_2,e_3,...} where e_i is a string
%! or a number. It is assumed that each element e_i is unique in set A_i.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item a
%! n*1 vector of doubles.
%! @item a
%! m*1 vector of doubles.
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item cprod
%! nm*2 matrix of doubles, the nodes cartesian product of sets @var{a} and @var{b}.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
%! @sp 2
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

[ F{1:nargin} ] = ndgrid( varargin{:} );

for i=1:nargin
    cprod(:,i) = F{i}(:);
end