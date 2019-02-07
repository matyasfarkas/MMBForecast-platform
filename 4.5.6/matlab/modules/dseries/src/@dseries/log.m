function ts = log(ts)

%@info:
%! @deftypefn {Function File} {@var{ts} =} log(@var{ts})
%! @anchor{log}
%! Apply the logarithm function to a Dynare time series object.
%!
%! @strong{Inputs}
%! @table @var
%! @item ts
%! Dynare time series object, instantiated by @ref{dseries}
%! @end table
%!
%! @strong{Outputs}
%! @table @var
%! @item ts
%! Dynare time series object with transformed data field.
%! @end table
%!
%! @strong{This function is called by:}
%! None.
%!
%! @strong{This function calls:}
%! None.
%!
%! @end deftypefn
%@eod:

% Copyright (C) 2011-2013 Dynare Team
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

if any(ts.data<eps)
    error('dseries::log: Input argument has to be strictly positive!')
end

for i=1:vobs(ts)
    ts.name(i) = {['log(' ts.name{i} ')']};
    ts.tex(i) = {['\log(' ts.tex{i} ')']};
end

ts.data = log(ts.data);