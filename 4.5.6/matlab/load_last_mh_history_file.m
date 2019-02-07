function info = load_last_mh_history_file(MetropolisFolder, ModelName)
% function info = load_last_mh_history_file(MetropolisFolder, ModelName)
% Loads the last mh_history_file
% Inputs:
%   MetropolisFolder    [char]      Name of the metropolis subfolder
%   ModelName           [char]      Name of the mod-file
% Outputs:
%   info                [struct]    structure storing the MH history
%
% Notes: The record structure is written to the caller workspace via an
% assignin statement.

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

% Consistency with older versions of Dynare.
if isequal(length(mh_history_files),0)
    if exist([BaseName '_mh_history.mat'])
        format_mh_history_file = 1; % old Dynare format
    else
        error(['Estimation::load_mh_file: I cannot find any mh-history file in ' MetropolisFolder '!'])
    end
else
    format_mh_history_file = 0;
end

if format_mh_history_file %needed to preserve backward compatibility
    load([BaseName '_mh_history.mat']);
    record.LastLogPost = record.LastLogLiK;
    record.InitialLogPost = record.InitialLogLiK;
    record.LastSeeds = record.Seeds;
    record.AcceptanceRatio = record.AcceptationRates;
    record.InitialSeeds = NaN; % This information is forever lost...
    record.MCMCConcludedSuccessfully = NaN; % This information is forever lost...
    record.MAX_nruns=NaN(size(record.MhDraws,1),1); % This information is forever lost...
    record = rmfield(record,'LastLogLiK');
    record = rmfield(record,'InitialLogLiK');
    record = rmfield(record,'Seeds');
    record = rmfield(record,'AcceptationRates');
    save([BaseName '_mh_history_0.mat'],'record');
else
    load([BaseName '_mh_history_' num2str(length(mh_history_files)-1) '.mat']);
    % add fields that have later been introduced
    if ~isfield(record,'MCMCConcludedSuccessfully')
        record.MCMCConcludedSuccessfully = NaN; % This information is forever lost...
    end
    if ~isfield(record,'MAX_nruns')
        record.MAX_nruns=NaN(size(record.MhDraws,1),1); % This information is forever lost...
    end
end

if isequal(nargout,0)
    assignin('caller', 'record', record);
else
    info = record;
end