function b = isolder(f, F) % --*-- Unitary tests --*--

% Returns true if f is older than any file in folder (and subfolders) F.
%
% INPUTS
% - f   [string]  file name
% - F   [string]  folder name
%
% OUTPUT
% - b   [logical]

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

if nargin<2 || isempty(F)
    F = pwd();
end

b = true;

files = dir(F);
tfile = dir(f);
tdate = tfile.datenum;

for i=1:length(files)
    if isequal(files(i).name, '.') || isequal(files(i).name, '..')
        continue
    end
    if isequal(files(i).name, f)
        continue
    end
    if files(i).isdir
        b = isolder(f, files(i).name);
        if ~b
            return
        end
        continue
    end
    b = tdate<files(i).datenum;
    if ~b
        return
    end
end

%@test:1
%$ t = false(3,1);
%$ mkdir toto;
%$ cd toto
%$ mkdir titi;
%$ mkdir tata;
%$ cd tata
%$ mkdir tutu;
%$ cd tutu
%$ system('touch a.m');
%$ cd ..
%$ system('touch b.m');
%$ system('touch c.m');
%$ cd ../titi
%$ system('touch d.m');
%$ cd ..
%$ pause(1)
%$ system('touch e.m');
%$ t(1) = isequal(isolder('e.m'), false);
%$ pause(1)
%$ system('touch tata/tutu/a.m');
%$ system('touch tata/b.m');
%$ system('touch tata/c.m');
%$ system('touch titi/d.m');
%$ t(2) = isequal(isolder('e.m'), true);
%$ pause(1)
%$ system('touch e.m');
%$ t(3) = isequal(isolder('e.m'), false);
%$ cd ..
%$ if isoctave()
%$   confirm_recursive_rmdir(false, 'local');
%$ end
%$ rmdir('toto','s');
%$ T = all(t);
%@eof:1
