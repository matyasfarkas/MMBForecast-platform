function A = subsref(A, S)
%function A = subsref(A, S)

% Copyright (C) 2013-2017 Dynare Team
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

switch S(1).type
  case '.'
    switch S(1).subs
      case fieldnames(A)
        A = A.(S(1).subs);
      case methods(A)
            if areParensNext(S)
                A = feval(S(1).subs, A, S(2).subs{:});
                S = shiftS(S,1);
            else
                A = feval(S(1).subs, A);
            end
      otherwise
        error(['@section.subsref: unknown field or method: ' S(1).subs]);
    end
  case '()'
    A = A.elements{S(1).subs{:}};
  case '{}'
    error(['@section.subsref: ' S(1).type ' indexing not supported.']);
  otherwise
    error('@section.subsref: impossible case')
end

S = shiftS(S,1);
if length(S) >= 1
    A = subsref(A, S);
end
end
