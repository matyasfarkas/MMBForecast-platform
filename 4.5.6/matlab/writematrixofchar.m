function str = writematrixofchar(m)

% Writes a matrix of char in a string.
%
% INPUTS
% - m   [char] matrix of char.
%
% OUTPUTS
% - str [char]
%
% EXAMPLE
% >> writematrixofchar(['a'; 'b'])
%
% ans =
%
% ['a'; 'b']
%
% where the returned argument is a string which can be evaluated or printed.

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

if ~(ischar(m) && ismatrix(m))
    error('Input has to be a matrix of char!')
end

str = '[';
for i=1:size(m, 1)
    str = sprintf('%s''%s''', str, m(i,:));
    if i<size(m, 1)
        str = sprintf('%s; ', str);
    end
end
str = sprintf('%s]', str);