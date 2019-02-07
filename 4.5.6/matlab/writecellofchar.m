function str = writecellofchar(c)

% Writes a two dimensional cell of char in a string.
%
% INPUTS
% - c   [cell] cell of char.
%
% OUTPUTS
% - str [char]
%
% EXAMPLES
% >> writecellofchar({'a', {'b'; 'c'}})
%
% ans =
%
% {'a', {'b'; 'c'}}
%
% >> writecellofchar({'a', ['b'; 'c'], 'd'})
%
% ans =
%
%{'a', '['b'; 'c']', 'd'}

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

str = '{';
for i=1:size(c, 1)
for j=1:size(c, 2)
if iscell(c{i,j})
str = sprintf('%s%s', str, writecellofchar(c{i, j}));
elseif ischar(c{i, j})
if size(c{i, j}, 1)>1
str = sprintf('%s''%s''', str, writematrixofchar(c{i, j}));
else
str = sprintf('%s''%s''', str, c{i, j});
end
else
error('Type not implemenented!')
end
if j<size(c, 2)
str = sprintf('%s, ', str);
end
end
if i<size(c, 1)
str = sprintf('%s; ', str);
end
end
str = sprintf('%s}', str);