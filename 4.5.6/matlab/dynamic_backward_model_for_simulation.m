function [r, J] = dynamic_backward_model_for_simulation(z, dynamicmodel, ylag, x, params, steady_state, it_)

% Copyright (C) 2017 Dynare Team
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

% Get indices of the variables appearing at time t.
% NOTE: It is assumed that all variables appear at time t in the model.
idy = length(ylag)+(1:length(z));

% Build y vector to be passed to the dynamic model.
y = zeros(length(ylag)+length(z), 1);
y(1:length(ylag)) = ylag;
y(idy) = z;

if nargout>1
    % Compute residuals and jacobian of the full dynamic model.
    [r, Jacobian] = feval(dynamicmodel, y, x, params, steady_state, it_);
else
    % Compute residuals and return.
    r = feval(dynamicmodel, y, x, params, steady_state, it_);
    return
end

% If the jacobian is computed, remove the columns related to the innovations
% and the variables appearing at time t-1.
J = Jacobian(:,idy);