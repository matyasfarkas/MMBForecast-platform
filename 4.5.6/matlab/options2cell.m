function c = options2cell(o)

% Converts an option structure as a cell of NAME and VALUE pairs.
%
% INPUTS
%  o o       matlab's structure holding a set of options (each field name is the name of an option and the associated content is the value of the option).
%
% OUTPUTS
%  o c       matlab's cell row array of the form {NAME1, VALUE1, NAME2, VALUE2, NAME3, VALUE3, ...}.

% Copyright (C) 2013-2017 Dynare Team.
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

s = fieldnames(o);
c = {};
j = 1;

for i=1:length(s)
    c(j) = {s{i}};
    c(j+1) = {o.(s{i})};
    j = j+2;
end