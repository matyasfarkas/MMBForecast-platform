function varargout = size(o, varargin)

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

switch nargout
  case 0
    size(o.data, varargin{:})
  case 1
    varargout{1} = size(o.data, varargin{:});
  case 2
    if isequal(nargin, 1)
        varargout{1} = size(o.data, 1);
        varargout{2} = size(o.data, 2);
    else
        error('dseries::size: Wrong calling sequence!')
    end
  otherwise
    error('dseries::size: Wrong calling sequence! Cannot return more than two arguments.')
end