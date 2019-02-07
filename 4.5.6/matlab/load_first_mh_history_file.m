function info = load_first_mh_history_file(MetropolisFolder, ModelName)

% This routine requires that the MCMC draws were obtained with a dynare version greater than 4.3.3.

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

% record is also a Matlab function.
record = 0;

% Get the list of all the mh_history files.
BaseName = [MetropolisFolder filesep ModelName];
mh_history_files = dir([BaseName '_mh_history_*.mat']);

if isequal(length(mh_history_files),0)
    error(['Estimation::load_mh_file: I cannot find any mh-history file in ' MetropolisFolder '!'])
end

load([BaseName '_mh_history_0.mat']);

if isequal(nargout,0)
    assignin('caller', 'record', record);
else
    info = record;
end