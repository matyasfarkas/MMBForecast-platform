function [shocks, spfm_exo_simul, innovations, DynareResults] = extended_path_shocks(innovations, ep, exogenousvariables, sample_size,DynareModel,DynareOptions, DynareResults)

% Copyright (C) 2016-2017 Dynare Team
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

% Simulate shocks.
if isempty(exogenousvariables)
    switch ep.innovation_distribution
      case 'gaussian'
        shocks = zeros(sample_size, DynareModel.exo_nbr);
        shocks(:,innovations.positive_var_indx) = transpose(transpose(innovations.covariance_matrix_upper_cholesky)*randn(innovations.effective_number_of_shocks,sample_size));
      case 'calibrated'
        options = DynareOptions;
        options.periods = options.ep.periods;
        oo = make_ex_(DynareModel,options,DynareResults);
        shocks = oo.exo_simul(2:end,:);
      otherwise
        error(['extended_path:: ' ep.innovation_distribution ' distribution for the structural innovations is not (yet) implemented!'])
    end
else
    shocks = exogenousvariables;
    innovations.positive_var_indx = find(sum(abs(shocks)>0));
end

% Copy the shocks in exo_simul
DynareResults.exo_simul = shocks;
spfm_exo_simul = repmat(DynareResults.exo_steady_state',ep.periods+2,1);