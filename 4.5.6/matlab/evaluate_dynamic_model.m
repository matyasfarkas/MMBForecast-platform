function residuals = evaluate_dynamic_model(dynamicmodel, endogenousvariables, exogenousvariables, params, steadystate, leadlagincidence, samplesize)

% Copyright (C) 2016 Dynare Team
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

ny = length(steadystate);
periods = rows(exogenousvariables);

residuals = zeros(ny,samplesize);
icols = find(leadlagincidence');

for t = 2:samplesize+1
    residuals(:,t-1) = dynamicmodel(endogenousvariables(icols), exogenousvariables, params, steadystate, t);
    icols = icols + ny;
end