function [DirectoryName, info] = CheckPath(type,dname)
% Creates the subfolder "./M_.dname/type" if it does not exist yet.
%
% INPUTS
%    type   [string]    Name of the subfolder.
%    dname  [string]    Name of the directory
%
% OUTPUTS
%    DirectoryName      string, name of the directory (with path).
%    info               integer scalar, equal to 1 if the routine creates directory dname/type, zero otherwise.
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2005-2017 Dynare Team
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

info = 0;

DirectoryName = [ dname '/' type ];

if ~isdir(dname)
    % Make sure there isn't a file with the same name, see trac ticket #47
    if isfile(dname)
        delete(dname)
    end
    mkdir('.', dname);
end

if ~isdir(DirectoryName)
    % Make sure there isn't a file with the same name, see trac ticket #47
    if isfile(DirectoryName)
        delete(DirectoryName)
    end
    mkdir('.',DirectoryName);
    info = 1;
end
