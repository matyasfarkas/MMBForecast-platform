function check_dsge_var_model(Model, EstimatedParameters, BayesInfo)

% Check if the dsge model can be estimated with the DSGE-VAR approach.

% Copyright (C) 2013-2014 Dynare Team
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

if EstimatedParameters.nvn
    error('Estimation::DsgeVarLikelihood: Measurement errors are not allowed!')
end

if EstimatedParameters.ncn
    error('Estimation::DsgeVarLikelihood: Measurement errors are not allowed!')
end

if any(vec(Model.H))
    error('Estimation::DsgeVarLikelihood: Measurement errors are not allowed!')
end

if EstimatedParameters.ncx
    error('Estimation::DsgeVarLikelihood: Structural innovations cannot be correlated using Dynare''s interface! Introduce the correlations in the model block instead.')
end

if Model.exo_nbr>1 && any(vec(tril(Model.Sigma_e,-1)))
    error('Estimation::DsgeVarLikelihood: Structural innovations cannot be correlated using Dynare''s interface! Introduce the correlations in the model block instead.')
end

if isequal(BayesInfo.with_trend,1)
    error('Estimation::DsgeVarLikelihood: Linear trend is not yet implemented!')
end
