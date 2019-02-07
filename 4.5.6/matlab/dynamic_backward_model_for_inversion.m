function [r, J] = dynamic_backward_model_for_inversion(z, dynamicmodel, ylag, ycur, x, params, steady_state, it_, ModelInversion)

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

% Set up y
y = zeros(length(ylag)+ModelInversion.nyfree+ModelInversion.nyctrl,1);
y(1:length(ylag)) = ylag;

y(ModelInversion.y_constrained_id) = ycur;
if ModelInversion.nyfree
    y(ModelInversion.y_free_id) = z(1:ModelInversion.nyfree);
end

% Update x
x(it_, ModelInversion.x_free_id) = transpose(z(ModelInversion.nyfree+(1:ModelInversion.nxfree)));

if nargout>1
    [r, Jacobian] = feval(dynamicmodel, y, x, params, steady_state, it_);
else
    r = feval(dynamicmodel, y, x, params, steady_state, it_);
    return
end

J = Jacobian(:,ModelInversion.J_id);