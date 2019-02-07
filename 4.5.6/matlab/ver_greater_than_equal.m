function tf = ver_greater_than_equal(ver1, ver2)
%function tf = ver_greater_than_equal(ver1, ver2)
% ver1 >= ver2 ? 1 : 0;
%
% INPUTS
%    ver1    [string]    software version number
%    ver2    [string]    software version number
%
% OUTPUTS
%    tf      [bool]      true if ver1 > ver2
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2015 Dynare Team
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

tf = ver_greater_than(ver1, ver2) || strcmp(ver1, ver2);
end