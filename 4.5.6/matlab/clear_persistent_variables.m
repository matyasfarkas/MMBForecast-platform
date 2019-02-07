function clear_persistent_variables(folder, writelistofroutinestobecleared)

% Clear all the functions with persistent variables in directory folder (and subdirectories).

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


if nargin<2
    writelistofroutinestobecleared = false;
end

if nargin<1 || isempty(folder)
    folder = pwd();
end

DYNARE_FOLDER = strrep(which('dynare'),'dynare.m','');

if writelistofroutinestobecleared
    if ~exist('list_of_functions_to_be_cleared.m') || isolder(sprintf('%slist_of_functions_to_be_cleared.m', DYNARE_FOLDER), DYNARE_FOLDER)
        if isunix() || ismac()
            [status, output] = system(sprintf('grep -lr ^persistent %s', folder));
            list_of_files = strsplit(output);
            list_of_files(find(cellfun(@isempty, list_of_files))) = [];
        else
            [status, output] = system(sprintf('findstr /B/S/M persistent %s\\*', folder));
            list_of_files = strsplit(output);
            list_of_files(find(cellfun(@isempty, list_of_files))) = [];
            i = 1; mobius = true;
            while mobius
                if i>length(list_of_files)
                    break
                end
                if ismember(list_of_files{i},{'FINDSTR:', 'ignored', '//'})
                    list_of_files(i) = [];
                else
                    i = i + 1;
                end
            end
        end
        [paths, list_of_functions, extensions] = cellfun(@fileparts, list_of_files, 'UniformOutput',false);
        cellofchar2mfile(sprintf('%slist_of_functions_to_be_cleared.m', DYNARE_FOLDER), list_of_functions)
    end
    return
end

list_of_functions_to_be_cleared;
clear(list_of_functions{:});