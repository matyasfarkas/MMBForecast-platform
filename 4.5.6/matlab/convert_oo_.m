function oo_ = convert_oo_(M_, options_, oo_, from_ver, to_ver)
%function oo_ = convert_oo_(M_, options_, oo_, from_ver, to_ver)
% Converts oo_ from oo_.dynare_version to ver
%
% INPUTS
%    M_          [struct]    dynare model struct
%    options_    [struct]    dynare options struct
%    oo_         [struct]    dynare output struct
%    from_ver    [string]    original oo_ output version
%    to_ver      [string]    desired oo_ output version
%
% OUTPUTS
%    oo_         [struct]    dynare output struct
%
% SPECIAL REQUIREMENTS
%    none

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
% along = with Dynare.  If not, see <http://www.gnu.org/licenses/>.

check_valid_ver(to_ver);
check_valid_ver(from_ver);

MIN_VER = '4.4';
MAX_VER = '4.5';

if ver_less_than(to_ver, MIN_VER)
    error(['Can only convert as far back as Dynare ' MIN_VER
           '. All versions before have the same oo_ structure.']);
end

if ver_greater_than(to_ver, MAX_VER)
    error(['Can only convert up to Dynare ' MAX_VER]);
end

min_ver = strsplit(MIN_VER, {'.', '-'});
max_ver = strsplit(MAX_VER, {'.', '-'});
from_ver_split = strsplit(from_ver, {'.', '-'});
to_ver_split = strsplit(to_ver, {'.', '-'});

if length(from_ver_split) ~= 2 || length(to_ver_split) ~= 2
    error('The version numbers may only be of the form X.Y');
end

if to_ver_split{2} > from_ver_split{2}
    new_to_ver = [to_ver_split{1} '.' num2str(str2double(to_ver_split{2})-1)];
else
    new_to_ver = [to_ver_split{1} '.' num2str(str2double(to_ver_split{2})+1)];
end

if strcmp(from_ver, to_ver)
    return
end

if ver_greater_than(to_ver, from_ver)
    moving_up = 1;
else
    moving_up = -1;
end

oo_ = convert_oo_(M_, options_, oo_, from_ver, new_to_ver)

if abs(to_ver_split{2} - (from_ver_split{2} - moving_up)) > 1
    new_from_ver = [to_ver_split{1} '.' num2str(str2double(to_ver_split{2}) - moving_up)];
else
    new_from_ver = from_ver;
end

eval(['oo_ = convert_dyn_' strrep(new_from_ver, '.', '') '_to_' ...
      strrep(to_ver, '.', '') '(M_, options_, oo_);']);
end
