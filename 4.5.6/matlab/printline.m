function varargout = printline(n, s, fid)
% This function print a line formed by replicating a symbol s.
%
% INPUTS
%
%   n  [integer]    Length of the printed line
%   s  [char]       Symbol used to draw the line (+, -, =, ...)
%   f  [integer]    file id returned by fopen
%
% OUTPUTS
%   None

% Copyright (C) 2015-2017 Dynare Team
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

if nargin<3
    f = 1;
    if nargin<2
        s = '-';
        if ~nargin
            error('printline: First argument is mandatory!')
        end
    end
end

S = s;

for i=2:n
    S = sprintf('%s%s',S,s);
end

if nargout
    varargout(1) = { sprintf('%s',S) };
else
    fprintf(f,sprintf('%s\n',S))
end