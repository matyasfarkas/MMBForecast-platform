function info = delete_mh_history_files(MetropolisFolder, ModelName)

% Copyright (C) 2013-2017 Dynare Team
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

BaseName = [MetropolisFolder filesep ModelName];

% Get the list of all the mh_history files.
mh_history_files = dir([BaseName '_mh_history_*.mat']);

% Delete the mh-history files (old or new format) if any.
if isequal(length(mh_history_files),0)
    if exist([BaseName '_mh_history.mat'])
        delete([BaseName '_mh_history.mat'])
    end
else
    delete([BaseName '_mh_history_*.mat'])
end