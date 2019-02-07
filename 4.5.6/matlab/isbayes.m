function l = isbayes(estim_params_)

% Returns true iff bayesian priors over parameters are defined.

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

l = false;

if ~isstruct(estim_params_)
    return
end

if isempty(estim_params_)
    return
end

ptypes = {'param_vals', 5; 'var_exo', 5 ; 'var_endo', 5; 'corrx', 6; 'corrn', 6};

for i=1:size(ptypes, 1)
    if isfield(estim_params_, ptypes{i, 1})
        tmp = estim_params_.(ptypes{i, 1});
        if ~isempty(tmp)
            if any(tmp(:,ptypes{i, 2})>0)
                l = true;
                return
            end
        end
    end
end