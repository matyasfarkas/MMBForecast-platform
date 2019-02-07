function info = is_unitary_test_available(fun)

% Decides if unitary tests defined in a matlab routine (file) have to be run
% by checking the content of the first line.
%
% INPUTS
%  - fun  [string], name of the routine (with full relative path)
%
% OUTPUTS
%  - info [integer], scalar equal to 1 if unitary tests must be run, 0 otherwise.

% Copyright (C) 2013-2017 Dynare Team
%
% This file is part of Dynare (m-unit-tests module).
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare's m-unit-tests module is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
% more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

info = 0;

fid = fopen(fun,'r');
first_line = fgetl(fid);
fclose(fid);

if strfind(first_line,'% --*-- Unitary tests --*--')
    info = 1;
end