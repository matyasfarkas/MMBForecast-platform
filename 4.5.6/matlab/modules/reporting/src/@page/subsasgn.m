function B = subsasgn(A, S, V)
% function B = subsasgn(A, S, V)

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

B = A;
if length(S) > 1
    for i=1:(length(S)-1)
        B = subsref(B, S(i));
    end
    B = subsasgn(B, S(end), V);
    B = subsasgn(A, S(1:(end-1)), B);
    return
end

switch S.type
  case '()'
    index = S.subs{:};
    assert(isnumeric(index));
    B{index} = V;
  case '.'
    switch S.subs
      case fieldnames(A)
        B.(S.subs) = V;
      otherwise
        error(['@page.subsasgn: field ' S.subs 'does not exist']);
    end
  otherwise
    error('@page.subsasgn: syntax error');
end
end