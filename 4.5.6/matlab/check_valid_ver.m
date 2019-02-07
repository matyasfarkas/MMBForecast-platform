function check_valid_ver(ver)
%function check_valid_ver(ver)
% Checks that ver is valid
%
% INPUTS
%    ver    [string]    dynare version number
%
% OUTPUTS
%    none
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

test_ver = strsplit(ver, {'.', '-'});
errmsg = 'check_valid_ver: the desired version must be in the proper format';

assert((length(test_ver) == 2 || length(test_ver) == 3) ...
       && ~isempty(str2double(test_ver{1})) ...
       && ~isempty(str2double(test_ver{2})), ...
       errmsg);
if length(test_ver) == 3 && isnan(str2double(test_ver{3}))
    assert(strcmp(test_ver{3}, 'unstable'), errmsg);
end
end
