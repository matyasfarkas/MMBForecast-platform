function overallacceptanceratio = compute_overall_acceptance_ratio(MetropolisFolder, ModelName)

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
mh_history_files = dir([BaseName '_mh_history_*.mat']);

n = length(mh_history_files);

load([BaseName '_mh_history_' num2str(0)]);
TotalNumberOfDraws = record.MhDraws(end,1);
TotalNumberOfAcceptedProposals = record.AcceptanceRatio*record.MhDraws(end,1);

for i=2:n
    load([BaseName '_mh_history_' num2str(i-1)]);
    TotalNumberOfDraws = TotalNumberOfDraws + record.MhDraws(end,1);
    TotalNumberOfAcceptedProposals = TotalNumberOfAcceptedProposals + record.AcceptanceRatio*record.MhDraws(end,1);
end

overallacceptanceratio = TotalNumberOfAcceptedProposals/TotalNumberOfDraws;