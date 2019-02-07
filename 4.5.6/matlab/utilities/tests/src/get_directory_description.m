function flist = get_directory_description(basedir)

% Lists recursively all the *.m files in a directory.
%
% INPUTS
%  - basedir [string], the name of the directory to be inspected.
%
% OUTPUTS
%  - flist   [cell of strings], the files under basedir (and subfolders).

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

if ~nargin || isempty(basedir)
    % Current directory is the default value of basedir
    basedir = '.';
end

dd = dir(basedir);
flist = {};
file = 1;

for f=1:length(dd)
    if ~(isequal(dd(f).name,'.') || isequal(dd(f).name,'..'))
        if dd(f).isdir
            r = get_directory_description([ basedir filesep dd(f).name]);
            flist = { flist{:} r{:} };
        else
            % Filter out files without m extension.
            [make, my, ext] = fileparts(dd(f).name);
            if isequal(ext, '.m')
                flist{length(flist)+1} = [basedir filesep dd(f).name];
            end
        end
        file = file + 1;
    end
end